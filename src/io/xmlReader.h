/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2010 Jonas Latt, Jonas Fietz, Mathias Krause
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

/** \file
 * Input/Output in XML format -- header file.
 */
#ifndef XML_IO_H
#define XML_IO_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <typeinfo>

#include <tinyxml.h>

#include "io/ostreamManager.h"
#include "communication/mpiManager.h"
#include "xmlReaderOutput.h"

namespace olb {

namespace util {
  template <class T, unsigned DIM> class ADf;
}

class XMLreader {
  friend class olb::XMLreaderOutput;
public:
  /**
   * Constructs a new XMLreader from another XMLreader
   * \param pParent The new root node for the XMLreader
   */
  XMLreader( TiXmlNode* pParent, OutputChannel outputChannel = OutputChannel::ERRCHANNEL);
  /// Constructs a new XMLreader from a XML file fName
  XMLreader( const std::string& fName, OutputChannel outputChannel = OutputChannel::ERRCHANNEL);
  /// destructor
  ~XMLreader();
  /// Prints out the XML structure read in, mostly for debugging purposes
  //void print(int indent) const;
  /**
   * Read a value from the xml file
   * \param reference to return the value
   * \return returns the value
   */
  //bool read(bool& value, bool verbose = true) const;
  template <typename T> bool read(T& value, bool verboseOn = true, bool exitIfMissing=false) const;
  template <typename T,unsigned DIM> bool read(util::ADf<T,DIM>& value, bool verboseOn = true, bool exitIfMissing=false) const;
  template <typename T> bool read(std::vector<T>& value, bool verboseOn = true, bool exitIfMissing=false) const;
  template <typename T> T get(bool verboseOn = true, bool exitIfMissing=false) const;
  /// This wrapper function reads the given parameter from the "type_parameter" and "name_parameter_1" or "name_parameter_2" tag and prints a warning, if the parameter can not be read.
  /// The warning contains the default value, if available. Will exit(1) if exitIfMissing == true. The warning is not displayed, if showWarning == false.
  template<typename ParameterType>
  bool readOrWarn(std::string name_parameter_1,
                   std::string name_parameter_2, std::string name_parameter_3,
                   ParameterType& var, bool defaultAvailable = true, bool exitIfMissing = false, bool showWarning = true) const;
  /// \return a Subtree placed at name \param name The name from which to take the subtree
  template<typename ParameterType>
  bool readOrWarn(std::string name_parameter_1,
                   std::string name_parameter_2, std::string name_parameter_3, std::string name_parameter_4,
                   ParameterType& var, bool defaultAvailable = true, bool exitIfMissing = false, bool showWarning = true) const;
  /// \return a Subtree placed at name \param name The name from which to take the subtree
  XMLreader const& operator[] (std::string name) const;
  /**
   * Returns an iterator.begin() of the child XMLreader
   * This means an iterator to the next level on an XML tree.
   */
  std::vector<XMLreader*>::const_iterator begin() const;
  /**
   * Returns an iterator.end() of the child XMLreader
   * This means an iterator to the next level on an XML tree.
   */
  std::vector<XMLreader*>::const_iterator end() const;
  /// switch warnings on/off
  void setWarningsOn(bool warnings) const;
  /// return the name of the element
  std::string getName() const;
  /// return the text of the element
  std::string getText() const;
  /// \return the value of attribute
  std::string getAttribute(const std::string& aName) const;

  /// handling all the output for the XMLreader
  XMLreaderOutput _output;
private:
  void mainProcessorIni(TiXmlNode* pParent);
  void slaveProcessorIni();
  XMLreader();
private:
  mutable bool _warningsOn;
  std::string _text;
  std::string _name;
  static XMLreader _notFound;
  OutputChannel _outputChannel;
protected:
  std::map<std::string, std::string> _attributes;
  std::vector<XMLreader*> _children;
};

// methods with template

template <typename T, unsigned DIM>
bool XMLreader::read(util::ADf<T,DIM>& value, bool verboseOn, bool exitIfMissing) const
{
  std::stringstream valueStr(_text);
  T tmp = T();
  if (!(valueStr >> tmp)) {
//    if ( _verboseOn ) {
//      clout << std::string("Error: cannot read value from XML element ") << _name << std::endl;
//    }
    _output.printWarning(_name, "ADf vector", "", verboseOn, exitIfMissing);
    return false;
  }
  value = util::ADf<T,DIM>(tmp);
  return true;
}

template <typename T>
bool XMLreader::read(std::vector<T>& values, bool verboseOn, bool exitIfMissing ) const
{
  std::stringstream multiValueStr(_text);
  std::string word;
  std::vector<T> tmp(values);
  while (multiValueStr>>word) {
    std::stringstream valueStr(word);
    T value;
    if (!(valueStr >> value)) {
//      if ( verboseOn ) {
//        clout << std::string("Error: cannot read value array from XML element ") << _name << std::endl;
//      }
      _output.printWarning(_name, "std::vector", "", verboseOn, exitIfMissing);
      return false;
    }
    tmp.push_back(value);
  }
  values.swap(tmp);
  return true;
}

template <typename T>
T XMLreader::get(bool verboseOn, bool exitIfMissing) const
{
  std::stringstream valueStr(_text);
  T tmp = T();
  if (!(valueStr >> tmp)) {
//    if ( verboseOn ) {
//      clout << "Error: cannot read value from XML element " << _name << std::endl;
//    }
    _output.printWarning(_name, typeid(T).name(), "", verboseOn, exitIfMissing);
  }
  return tmp;
}


template<typename ParameterType>
bool XMLreader::readOrWarn(std::string name_parameter_1,
                   std::string name_parameter_2, std::string name_parameter_3,
                   ParameterType& var, bool defaultAvailable, bool exitIfMissing, bool showWarning) const
{
  // deactivate default warnings and show default values instead
  setWarningsOn(false);
  if (name_parameter_3 == "") {
    if (!(*this)[name_parameter_1][name_parameter_2].read<ParameterType>(var, false)) {
      _output.parameterReading({name_parameter_1, name_parameter_2}, var, defaultAvailable, exitIfMissing, showWarning);
      return false;
    }
    return true;
  }
  else{
    if (!(*this)[name_parameter_1][name_parameter_2][name_parameter_3].read<ParameterType>(var, false)) {
      _output.parameterReading({name_parameter_1, name_parameter_2, name_parameter_3}, var, defaultAvailable, exitIfMissing, showWarning);
      return false;
    }
    return true;
  }
  // turn default warnings on again
  setWarningsOn(true);
}

template<typename ParameterType>
bool XMLreader::readOrWarn(std::string name_parameter_1,
                   std::string name_parameter_2, std::string name_parameter_3, std::string name_parameter_4,
                   ParameterType& var, bool defaultAvailable, bool exitIfMissing, bool showWarning) const
{
  // deactivate default warnings and show default values instead
  setWarningsOn(false);
  if (name_parameter_3 == "") {
    if (!(*this)[name_parameter_1][name_parameter_2].read<ParameterType>(var, false)) {
      _output.parameterReading({name_parameter_1, name_parameter_2}, var, defaultAvailable, exitIfMissing, showWarning);
      return false;
    }
    return true;
  }
  else if(name_parameter_4 == ""){
    if (!(*this)[name_parameter_1][name_parameter_2][name_parameter_3].read<ParameterType>(var, false)) {
      _output.parameterReading({name_parameter_1, name_parameter_2, name_parameter_3}, var, defaultAvailable, exitIfMissing, showWarning);
      return false;
    }
    return true;
  }
  else {
    if (!(*this)[name_parameter_1][name_parameter_2][name_parameter_3][name_parameter_4].read<ParameterType>(var, false)) {
      _output.parameterReading({name_parameter_1, name_parameter_2,name_parameter_3,name_parameter_4},
                                    var, defaultAvailable, exitIfMissing, showWarning);
      return false;
    }
    return true;
  }
  // turn default warnings on again
  setWarningsOn(true);
}

XMLreader XMLreader::_notFound;

XMLreader::XMLreader()
{
  _name = "XML node not found";
  _warningsOn = true;
  _output = XMLreaderOutput(OutputChannel::ERRCHANNEL);
}

XMLreader::XMLreader( TiXmlNode* pParent, OutputChannel outputChannel) : _output(outputChannel)
{
  _outputChannel = outputChannel;
  _warningsOn = true;

  if (singleton::mpi().isMainProcessor()) {
    mainProcessorIni(pParent);
  }
  else {
    slaveProcessorIni();
  }
}

XMLreader::XMLreader(const std::string& fName, OutputChannel outputChannel) : _output(outputChannel)
{
  _outputChannel = outputChannel;
  _warningsOn = true;

  TiXmlDocument* doc = nullptr;
  int loadOK = false;
#ifdef PARALLEL_MODE_MPI  // parallel program execution
  if (singleton::mpi().isMainProcessor()) {
#endif
    std::string docName = std::string(fName);  // call copy constructor
    doc = new TiXmlDocument(docName.c_str());
    loadOK = doc->LoadFile();
    _output.loadFile(loadOK, fName);
#ifdef PARALLEL_MODE_MPI  // parallel program execution
  }
  if (singleton::mpi().isMainProcessor()) {
#endif
    mainProcessorIni(doc);
    delete doc;
#ifdef PARALLEL_MODE_MPI  // parallel program execution
  }
  else {
    slaveProcessorIni();
  }
#endif
}

XMLreader::~XMLreader()
{
  for (unsigned int iNode=0; iNode<_children.size(); ++iNode) {
    delete _children[iNode];
  }
}

void XMLreader::mainProcessorIni( TiXmlNode* pParent )
{
  assert (pParent->Type()==TiXmlNode::TINYXML_DOCUMENT || pParent->Type()==TiXmlNode::TINYXML_ELEMENT );
  if (pParent->Type() == TiXmlNode::TINYXML_DOCUMENT) {
    // ignore the surrounding PARAM-block
    pParent = pParent->FirstChildElement();
  }

  _name = pParent->ValueStr();
#ifdef PARALLEL_MODE_MPI  // parallel program execution
  singleton::mpi().bCast(&_name,1);
#endif

  TiXmlAttribute* attr = pParent->ToElement()->FirstAttribute();
  while (attr != nullptr) {
#ifdef PARALLEL_MODE_MPI  // parallel program execution
    int size = 0;
    std::string* key = const_cast<std::string*>(&attr->NameTStr());
    singleton::mpi().bCast(key, size);
    std::string* value = const_cast<std::string*>(&attr->ValueStr());
    singleton::mpi().bCast(value, size);
#endif
    _attributes[attr->NameTStr()] = attr->ValueStr();
    attr = attr->Next();
  }
#ifdef PARALLEL_MODE_MPI  // parallel program execution
  std::string tmpstr = "";
  int size = 0;
  singleton::mpi().bCast(&tmpstr, size);
  singleton::mpi().bCast(&tmpstr, size);
#endif


  TiXmlNode * pChild;
  int type = 0;
  for ( pChild = pParent->FirstChild(); pChild != nullptr; pChild = pChild->NextSibling()) {
    type = pChild->Type();
#ifdef PARALLEL_MODE_MPI  // parallel program execution
    singleton::mpi().bCast(&type, 1);
#endif
    if ( type==TiXmlNode::TINYXML_ELEMENT ) {
      _children.push_back( new XMLreader( pChild , _outputChannel) );
    }
    else if ( type==TiXmlNode::TINYXML_TEXT ) {
      _text = pChild->ToText()->ValueStr();
#ifdef PARALLEL_MODE_MPI  // parallel program execution
      singleton::mpi().bCast(&_text,1);
#endif
    }
  }
  type = TiXmlNode::TINYXML_UNKNOWN;
#ifdef PARALLEL_MODE_MPI  // parallel program execution
  singleton::mpi().bCast(&type, 1);
#endif
}

void XMLreader::slaveProcessorIni()
{
#ifdef PARALLEL_MODE_MPI  // parallel program execution

  singleton::mpi().bCast(&_name,1);
  std::string key = "";
  std::string value = "";
  int size = int();
  do {
    singleton::mpi().bCast(&key, size);
    singleton::mpi().bCast(&value, size);
    _attributes[key] = value;
  }
  while (key != "");
#endif

  int type=0;
  do {
#ifdef PARALLEL_MODE_MPI  // parallel program execution
    singleton::mpi().bCast(&type, 1);
#endif
    if ( type==TiXmlNode::TINYXML_ELEMENT ) {
      _children.push_back( new XMLreader( nullptr, _outputChannel ) );
    }
    else if ( type==TiXmlNode::TINYXML_TEXT ) {
#ifdef PARALLEL_MODE_MPI  // parallel program execution
      singleton::mpi().bCast(&_text,1);
#endif
    }
  }
  while (type != TiXmlNode::TINYXML_UNKNOWN);
}

XMLreader const& XMLreader::operator[] (std::string fName) const
{
  for (unsigned int iNode=0; iNode<_children.size(); ++iNode) {
    if (_children[iNode]->_name == fName) {
      return *_children[iNode];
    }
  }
  _output.readValue(_warningsOn, _name, fName);
  return _notFound;
}

std::vector<XMLreader*>::const_iterator XMLreader::begin() const
{
  return _children.begin();
}

std::vector<XMLreader*>::const_iterator XMLreader::end() const
{
  return _children.end();
}

std::string XMLreader::getName() const
{
  return _name;
}

std::string XMLreader::getText() const
{
  return _text;
}

void XMLreader::setWarningsOn(bool warnings) const
{
  _warningsOn = warnings;
  for (unsigned int iNode=0; iNode<_children.size(); ++iNode) {
    _children[iNode]->setWarningsOn(warnings);
  }
}

// template specialization for T=bool
template <>
bool XMLreader::read<bool>(bool& value, bool verboseOn, bool exitIfMissing) const
{
  std::stringstream valueStr(_text);
  std::string word;
  valueStr >> word;
  // Transform to lower-case, so that "true" and "false" are case-insensitive.
  std::transform(word.begin(), word.end(), word.begin(), ::tolower);
  if (!word.compare("true") || (word=="1")) {
    value = true;
    return true;
  }
  else if (!word.compare("false") || (word=="0")) {
    value=false;
    return true;
  }
  else {
    if ( verboseOn ) {
      std::stringstream ss;
      ss << ( value ? "true" : "false" );
      _output.printWarning(_name, "bool", ss.str(),  verboseOn, exitIfMissing);
    }
  }
  return false;
}

// template specialization for T=int
template <>
bool XMLreader::read<int>(int& value, bool verboseOn, bool exitIfMissing) const
{
  std::stringstream valueStr(_text);
  int tmp = int();
  if (!(valueStr >> tmp)) {
    std::stringstream ss;
    ss << value;
    _output.printWarning(_name, "int", ss.str(), verboseOn, exitIfMissing);
    return false;
  }
  value = tmp;
  return true;
}

// template specialization for T=std::size_t
template <>
bool XMLreader::read<std::size_t>(std::size_t& value, bool verboseOn, bool exitIfMissing) const
{
  std::stringstream valueStr(_text);
  std::size_t tmp = std::size_t();
  if (!(valueStr >> tmp)) {
    std::stringstream ss;
    ss << value;
    _output.printWarning(_name, "std::size_t", ss.str(), verboseOn, exitIfMissing);
    return false;
  }
  value = tmp;
  return true;
}

// template specialization for T=double
template <>
bool XMLreader::read<double>(double& value, bool verboseOn, bool exitIfMissing) const
{
  std::stringstream valueStr(_text);
  double tmp = double();
  if (!(valueStr >> tmp)) {
    _output.printWarning(_name, "double", std::to_string(value), verboseOn, exitIfMissing);
    return false;
  }
  value = tmp;
  return true;
}

// template specialization for T=long double
template <>
bool XMLreader::read<long double>(long double& value, bool verboseOn, bool exitIfMissing) const
{
  std::stringstream valueStr(_text);
  std::string tmp {};
  if (!(valueStr >> tmp)) {
    _output.printWarning(_name, "long double", std::to_string(value), verboseOn, exitIfMissing);
    return false;
  }
  value = std::stold(tmp);
  return true;
}

template <>
bool XMLreader::read<float>(float& value, bool verboseOn, bool exitIfMissing) const
{
  std::stringstream valueStr(_text);
  float tmp = float();
  if (!(valueStr >> tmp)) {
    std::stringstream ss;
    ss << value;
    _output.printWarning(_name, "float", ss.str(), verboseOn, exitIfMissing);
    return false;
  }
  value = tmp;
  return true;
}

template <>
bool XMLreader::read<std::string>(std::string& entry, bool verboseOn, bool exitIfMissing) const
{
  if (_name == "XML node not found") {
    return false;
  }
  std::stringstream valueStr(_text);
  std::string tmp = std::string();
  if (!(valueStr >> tmp)) {
    std::stringstream ss;
    ss << entry;
    _output.printWarning(_name, "string", ss.str(), verboseOn, exitIfMissing);
    return false;
  }

  entry = _text;
  return true;
}

std::string XMLreader::getAttribute(const std::string& aName) const
{
  std::map<std::string, std::string>::const_iterator it = _attributes.find(aName);
  if ( it == _attributes.end()) {
    return "Attribute not found.";
  }
  return it->second;
  //return attributes[aName];
}


}  // namespace olb

#endif  // XML_IO_H
