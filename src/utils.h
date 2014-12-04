#ifndef CAP_UTILS_H
#define CAP_UTILS_H

#include <boost/property_tree/ptree.hpp>
#include <string>
#include <vector>

namespace cap {

template <typename T>
std::vector<T> to_vector(std::string const & s);

template <typename T>
std::string to_string(std::vector<T> const & v);

//template <typename T>
//void traverse_recursive
//    ( boost::property_tree::ptree const &            parent
//    , boost::property_tree::ptree::path_type const & child_path
//    , boost::property_tree::ptree const &            child
//    , T & method
//    )
//{
//    method(parent, child_path, child);
//    for (boost::property_tree::ptree::const_iterator it = child.begin(); it != child.end(); ++it) {
//        boost::property_tree::ptree::path_type current_path = child_path / boost::property_tree::ptree::path_type(it->first);
//        traverse_recursive(parent, current_path, it->second, method);
//    }
//}
//
//template <typename T>
//void traverse(const boost::property_tree::ptree & parent, T & method)
//{
//  traverse_recursive(parent, "", parent, method);
//}
//
//void merge
//    ( boost::property_tree::ptree /* const */ &            parent
//    , boost::property_tree::ptree::path_type const & child_path
//    , boost::property_tree::ptree const &            child
//    ) 
//{
//    parent.put(child_path, child.data());
//}   

} // end namespace cap

#endif // CAP_UTILS_H
