# Check compiler version
SET(INTEL_COMPILER 1)
IF ( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0 )
MESSAGE(FATAL_ERROR "Require intel 8.0 or higher ")
ENDIF()

# Enable OpenMP
SET(ENABLE_OPENMP 1)
IF ( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 16 )
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -openmp")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
ELSE()
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -qopenmp")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
ENDIF()

# Set the std
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -g -debug inline-debug-info -std=c99")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -debug inline-debug-info")

# Set intel specfic flags (which we always want)
ADD_DEFINITIONS( -DADD_ )
ADD_DEFINITIONS( -DINLINE_ALL=inline )

# Suppress compile warnings
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -Wno-deprecated")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")

# Set extra optimization specific flags
SET( CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -restrict -unroll -ip" )
SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -restrict -unroll -ip" )
SET( CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}   -restrict -unroll -ip" )
SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -restrict -unroll -ip" )
SET( CMAKE_C_FLAGS_RELWITHDEBINFO   "${CMAKE_C_FLAGS_RELWITHDEBINFO}   -restrict -unroll -ip" )
SET( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -restrict -unroll -ip" )

# Use deprecated options prior to 11.1
SET(ICC_DEPRECATED_OPTS FALSE)
IF ( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11.1 )
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -prefetch ")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -prefetch" )
ELSEIF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 16 )  
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -opt-prefetch" )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -opt-prefetch" )
ELSE()
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -qopt-prefetch" )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-prefetch" )
ENDIF()

#check if -ftz is accepted
CHECK_C_COMPILER_FLAG( "${CMAKE_CXX_FLAGS} -ftz" INTEL_FTZ )
IF( INTEL_FTZ)
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -ftz" )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftz" )
ENDIF( INTEL_FTZ)

exec_program(grep  
ARGS flags /proc/cpuinfo
OUTPUT_VARIABLE CPU_FLAGS
RETURN_VALUE CPUINFO
)

# SSE4.2 option is available for 11.1  and higher
# not all the flags are used by the code
SET(HAVE_SSE42 0)
SET(HAVE_SSE41 0)
SET(HAVE_SSSE3 0)
SET(HAVE_SSE3 0)
SET(HAVE_SSE2 0)
SET(HAVE_SSE 0)
SET(SET_HAVE_SSE_FLAGS 0 CACHE BOOL "Use einspline routines with SSE intrinsics")

#check if the user has already specified -x option for cross-compiling.
if(CMAKE_CXX_FLAGS MATCHES "-x" OR CMAKE_C_FLAGS MATCHES "-x" OR
    CMAKE_CXX_FLAGS MATCHES "-ax" OR CMAKE_C_FLAGS MATCHES "-ax")
  # make sure that the user specifies -x for both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS.
  if(CMAKE_CXX_FLAGS MATCHES "-x" AND CMAKE_C_FLAGS MATCHES "-x")
    SET(SSE_OPT_SET TRUE)
  else() #(CMAKE_CXX_FLAGS MATCHES "-x" AND CMAKE_C_FLAGS MATCHES "-x")
    if(CMAKE_CXX_FLAGS MATCHES "-ax" AND CMAKE_C_FLAGS MATCHES "-ax")
      SET(SSE_OPT_SET TRUE)
    else()
      MESSAGE(FATAL_ERROR "if -xcode is specified by the user, it should be added in both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS!")
    endif()
  endif() #(CMAKE_CXX_FLAGS MATCHES "-x" AND CMAKE_C_FLAGS MATCHES "-x")
else() #(CMAKE_CXX_FLAGS MATCHES "-x" OR CMAKE_C_FLAGS MATCHES "-x")
  #check if -xHost is accepted
  CHECK_C_COMPILER_FLAG( "-xHost" INTEL_CC_FLAGS )
  IF(INTEL_CC_FLAGS)
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -xHost")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost")
    SET(SSE_OPT_SET TRUE)
  ELSE(INTEL_CC_FLAGS)
    SET(SSE_OPT_SET FALSE)
  ENDIF(INTEL_CC_FLAGS)
endif() #(CMAKE_CXX_FLAGS MATCHES "-x" OR CMAKE_C_FLAGS MATCHES "-x")

#set all the HAVE_XXX flags based on /proc/cpuinfo
#set the compiler option -xcode based on /proc/cpuinfo if -xhost is not accepted.

if(CPU_FLAGS MATCHES "avx2")
IF(NOT SSE_OPT_SET)
  IF(ICC_DEPRECATED_OPTS)
    MESSAGE(WARNING "AVX2 needs version 14.0 and higher.")
  ELSE(ICC_DEPRECATED_OPTS)
    SET(SSE_OPT_SET TRUE)
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -xCORE-AVX2")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xCORE-AVX2")
  ENDIF(ICC_DEPRECATED_OPTS)
ENDIF(NOT SSE_OPT_SET)
endif(CPU_FLAGS MATCHES "avx2")

if(CPU_FLAGS MATCHES "avx")
IF(NOT SSE_OPT_SET)
  IF(ICC_DEPRECATED_OPTS)
    MESSAGE(WARNING "AVX needs version 11.1 and higher.")
  ELSE(ICC_DEPRECATED_OPTS)
    SET(SSE_OPT_SET TRUE)
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -xAVX")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xAVX")
  ENDIF(ICC_DEPRECATED_OPTS)
ENDIF(NOT SSE_OPT_SET)
endif(CPU_FLAGS MATCHES "avx")

if(CPU_FLAGS MATCHES "sse4_2")
IF(SET_HAVE_SSE_FLAGS)
  SET(HAVE_SSE42 1)
ENDIF(SET_HAVE_SSE_FLAGS)
IF(NOT SSE_OPT_SET)
  IF(ICC_DEPRECATED_OPTS)
    MESSAGE(WARNING "SSE4.2 needs version 11.1 and higher.")
  ELSE(ICC_DEPRECATED_OPTS)
    SET(SSE_OPT_SET TRUE)
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -xSSE4.2")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xSSE4.2")
  ENDIF(ICC_DEPRECATED_OPTS)
ENDIF(NOT SSE_OPT_SET)
endif(CPU_FLAGS MATCHES "sse4_2")

if(CPU_FLAGS MATCHES "sse4_1")
IF(SET_HAVE_SSE_FLAGS)
  SET(HAVE_SSE41 1)
ENDIF(SET_HAVE_SSE_FLAGS)
IF(NOT SSE_OPT_SET)
  SET(SSE_OPT_SET TRUE)
  IF(ICC_DEPRECATED_OPTS)
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -xS")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xS")
  ELSE(ICC_DEPRECATED_OPTS)
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -xSSE4.1")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xSSE4.1")
  ENDIF(ICC_DEPRECATED_OPTS)
ENDIF(NOT SSE_OPT_SET)
endif(CPU_FLAGS MATCHES "sse4_1")

if(CPU_FLAGS MATCHES "ssse3")
IF(SET_HAVE_SSE_FLAGS)
  SET(HAVE_SSSE3 1)
ENDIF(SET_HAVE_SSE_FLAGS)
IF(NOT SSE_OPT_SET)
  SET(SSE_OPT_SET TRUE)
  IF(ICC_DEPRECATED_OPTS)
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -xT")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xT")
  ELSE(ICC_DEPRECATED_OPTS)
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -xSSSE3")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xSSSE3")
  ENDIF(ICC_DEPRECATED_OPTS)
ENDIF(NOT SSE_OPT_SET)
endif(CPU_FLAGS MATCHES "ssse3")

if(CPU_FLAGS MATCHES "sse3")
IF(SET_HAVE_SSE_FLAGS)
  SET(HAVE_SSE3 1)
ENDIF(SET_HAVE_SSE_FLAGS)
IF(NOT SSE_OPT_SET)
  SET(SSE_OPT_SET TRUE)
  IF(ICC_DEPRECATED_OPTS)
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -xP")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xP")
  ELSE(ICC_DEPRECATED_OPTS)
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -xSSE3")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xSSE3")
  ENDIF(ICC_DEPRECATED_OPTS)
ENDIF(NOT SSE_OPT_SET)
endif(CPU_FLAGS MATCHES "sse3")

if(CPU_FLAGS MATCHES "sse2")
IF(SET_HAVE_SSE_FLAGS)
  SET(HAVE_SSE2 1)
  SET(HAVE_SSE 1)
ENDIF(SET_HAVE_SSE_FLAGS)
IF(NOT SSE_OPT_SET)
  SET(SSE_OPT_SET TRUE)
  IF(ICC_DEPRECATED_OPTS)
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -xN")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xN")
  ELSE(ICC_DEPRECATED_OPTS)
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -xSSE2")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xSSE2")
  ENDIF(ICC_DEPRECATED_OPTS)
ENDIF(NOT SSE_OPT_SET)
endif(CPU_FLAGS MATCHES "sse2")

######################################################
#KCC needs to be used to build static libraries
######################################################
#set(CMAKE_AR xild) 
#set(CMAKE_CXX_CREATE_STATIC_LIBRARY "<CMAKE_AR> -lib -o <TARGET> <OBJECTS>")

