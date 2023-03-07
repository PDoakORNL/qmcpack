//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

#include "InputSection.h"
#include "Message/UniformCommunicateError.h"
#include "ModernStringUtils.hpp"
#include "Utilities/string_utils.h"

namespace qmcplusplus
{

void InputSection::readAttributes(xmlNodePtr cur, bool consume_name)
{
  std::string ename{lowerCase(castXMLCharToChar(cur->name))};
  std::string par_name(lowerCase(getXMLAttributeValue(cur, "name")));    
  xmlAttrPtr att = cur->properties;
  while (att != NULL)
  {
    // unsafe att->name is an xmlChar, xmlChar is a UTF-8 byte
    std::string name{lowerCase(castXMLCharToChar(att->name))};
    // issue here is that we don't want to consume the name of the parameter as that has a special status in a parameter tag.
    // This is due to the <parameter name="parameter_name>  == <parametere_name> tag equivalence :(
    if (!consume_name && name == "name")
    {
      att = att->next;
      continue;
    }
    if (!(isAttribute(name) || isCustomAttribute(name)))
    {
      std::stringstream error;
      error << "InputSection::readXML name " << name << " is not an attribute of " << section_name << "\n";
      throw UniformCommunicateError(error.str());
    }
    std::istringstream stream(castXMLCharToChar(att->children->content));
    if (isCustomAttribute(name))
    {
      std::string ename{lowerCase(castXMLCharToChar(cur->name))};
      std::string par_name(lowerCase(getXMLAttributeValue(cur, "name")));
      if (ename == "parameter" && par_name.size() > 1)
        ename = par_name;
      setFromStreamCustom(ename, name, stream);
    }
    else
      setFromStream(name, stream);
    att = att->next;
  }
}

void InputSection::readXML(xmlNodePtr cur)
{
  // read attributes both the root node and parameters can have these.
  // But its limited what you can do with attributes of parameter nodes

  readAttributes(cur, true);

  // read parameters
  xmlNodePtr element = cur->xmlChildrenNode;
  while (element != NULL)
  {
    std::string ename{lowerCase(castXMLCharToChar(element->name))};
    std::string name(lowerCase(getXMLAttributeValue(element, "name")));

    // If the name attribute is empty and the element name is not parameter
    // set name to the tag name. It's contents will get this name as key.
    if (name.size() < 1 && ename != "parameter")
      name = ename;
    // If the element is a "parameter" then its name is what the parameter actually is
    // So set ename to name.
    if (ename == "parameter" && name.size() > 1)
      ename = name;

    if (isDelegate(ename))
    {
      assert(delegate_factories_.find(ename) != delegate_factories_.end());
      std::string value_key;
      std::any value = delegate_factories_[ename](element, value_key);
      if (has(value_key) && !isMultiple(value_key))
        throw UniformCommunicateError("Input is invalid  " + section_name + " contains " + ename +
                                      " node with duplicate name " + value_key + "!");
      if (isMultiple(value_key))
      {
        if (has(value_key))
        {
          auto* value_vector = std::any_cast<std::vector<std::any>>(&(values_[value_key]));
          value_vector->push_back(value);
        }
        else
        {
          values_[value_key] = std::vector<std::any>{value};
        }
      }
      else
      {
        values_[value_key] = value;
      }
    }
    else if (isCustomParameter(ename))
    {
      std::istringstream stream(XMLNodeString{element});
      setFromStreamCustom(ename, name, stream);
      readAttributes(element, false);
    }
    else if (isParameter(ename))
    {
      std::istringstream stream(XMLNodeString{element});
      setFromStream(name, stream);
      readAttributes(element, false);
    }
    else if (ename != "text")
    {
      std::stringstream error;
      error << "InputSection::readXML node name " << ename << " is not handled by " << section_name << "\n";
      throw UniformCommunicateError(error.str());
    }
    element = element->next;
  }

  // assign default values for optional variables
  setDefaults();

  // check input validity
  checkValid();
  //report();
}


void InputSection::init(const std::unordered_map<std::string, std::any>& init_values)
{
  // assign inputted values
  for (auto& [name, value] : init_values)
    setFromValue(name, value);

  // assign default values for optional variables
  setDefaults();

  // check input validity
  checkValid();
  //report();
}

void InputSection::registerDelegate(const std::string& tag,
                                    std::function<std::any(xmlNodePtr cur, std::string& value_key)> factory)
{
  delegate_factories_[tag] = factory;
}

void InputSection::setDefaults()
{
  for (auto& [name, default_value] : default_values)
    if (!has(name))
      setFromValue(name, default_value);
}

void InputSection::setFromStream(const std::string& name, std::istringstream& svalue)
{
  if (isString(name) || isEnumString(name))
  {
    std::string value;
    svalue >> value;
    values_[name] = value;
  }
  else if (isMultiString(name))
  {
    std::vector<std::string> string_values;
    for (std::string value; svalue >> value;)
      string_values.push_back(value);
    values_[name] = string_values;
  }
  else if (isMultiReal(name))
  {
    std::vector<Real> real_values;
    for (Real value; svalue >> value;)
      real_values.push_back(value);
    values_[name] = real_values;
  }
  else if (isBool(name))
  {
    std::string sval;
    svalue >> sval;
    bool value    = sval == "yes" || sval == "true" || sval == "1";
    values_[name] = value;
  }
  else if (isInteger(name))
  {
    int value;
    svalue >> value;
    values_[name] = value;
  }
  else if (isReal(name))
  {
    Real value;
    svalue >> value;
    values_[name] = value;
  }
  else if (isPosition(name))
  {
    Position value;
    svalue >> value;
    values_[name] = value;
  }
  else
  {
    std::stringstream error;
    error << "InputSection::set_from_stream name " << name << " in " << section_name << " does not have a type\n";
    throw UniformCommunicateError(error.str());
  }
}

template<typename T>
void InputSection::setFromValue(const std::string& name, const T& value)
{
  if (isString(name) || isEnumString(name))
    values_[name] = std::any_cast<std::string>(value);
  else if (isMultiString(name))
    values_[name] = (std::any_cast<std::vector<std::string>>(value));
  else if (isMultiString(name))
    values_[name] = (std::any_cast<std::vector<Real>>(value));
  else if (isBool(name))
    values_[name] = std::any_cast<bool>(value);
  else if (isInteger(name))
    values_[name] = std::any_cast<int>(value);
  else if (isReal(name))
    values_[name] = std::any_cast<Real>(value);
  else if (isPosition(name))
    values_[name] = std::any_cast<Position>(value);
  else
  {
    std::stringstream error;
    error << "InputSection::set_from_value name " << name << " in " << section_name << " does not have a type\n";
    throw UniformCommunicateError(error.str());
  }
}

void InputSection::checkValid()
{
  // check that all required inputs are present
  for (auto& name : required)
    if (!has(name))
    {
      std::stringstream error;
      error << "InputSection::check_valid required variable " << name << " in " << section_name
            << " has not been assigned\n";
      throw UniformCommunicateError(error.str());
    }
  checkParticularValidity();
};

void InputSection::report() const
{
  auto& out = app_log();
  report(out);
}

void InputSection::report(std::ostream& out) const
{
  out << "\n" << section_name;
  for (auto& [name, value] : values_)
  {
    out << "\n  " << name << " = ";
    if (isString(name))
      out << std::any_cast<std::string>(value);
    else if (isBool(name))
      out << std::any_cast<bool>(value);
    else if (isInteger(name))
      out << std::any_cast<int>(value);
    else if (isReal(name))
      out << std::any_cast<Real>(value);
  }
  out << "\n\n";
}


std::any InputSection::lookupAnyEnum(const std::string& enum_name,
                                     const std::string& enum_value,
                                     const std::unordered_map<std::string, std::any>& enum_map)
{
  std::string enum_value_str(lowerCase(enum_name + "-" + enum_value));
  try
  {
    return enum_map.at(enum_value_str);
  }
  catch (std::out_of_range& oor_exc)
  {
    std::throw_with_nested(UniformCommunicateError("bad_enum_tag_value: " + enum_value_str));
  }
}

} // namespace qmcplusplus
