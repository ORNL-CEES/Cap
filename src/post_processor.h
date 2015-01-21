#ifndef CAP_POSTPROCESSOR_H
#define CAP_POSTPROCESSOR_H

#include <cap/mp_values.h>
#include <cap/boundary_values.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/sparse_matrix.h>
//#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <boost/property_tree/ptree.hpp>
#include <memory>
#include <unordered_map>

namespace cap {

//////////////////////// POSTPROCESSOR PARAMETERS ////////////////////////////
template <int dim>
class PostprocessorParameters { 
public:
    PostprocessorParameters(std::shared_ptr<boost::property_tree::ptree const> d)
        : database(d)
    { }
    virtual ~PostprocessorParameters() { }

    dealii::DoFHandler<dim> const * dof_handler;
    dealii::BlockVector<double> const * solution;

    std::shared_ptr<MPValues<dim> const>       mp_values;
    std::shared_ptr<BoundaryValues<dim> const> boundary_values;

    std::shared_ptr<boost::property_tree::ptree const> database;
};

//////////////////////// POSTPROCESSOR /////////////////////
template <int dim>
class Postprocessor {
public:
    Postprocessor(std::shared_ptr<PostprocessorParameters<dim> const> parameters);
    virtual ~Postprocessor() { }
    virtual void reset(std::shared_ptr<PostprocessorParameters<dim> const> parameters) { }

    dealii::Vector<double> const & get(std::string const & key) const;
    void get(std::string const & key, double & value) const;
    std::vector<std::string> get_vector_keys() const;
    

protected:
    dealii::DoFHandler<dim> const & dof_handler;
    dealii::BlockVector<double> const & solution;

    std::shared_ptr<MPValues<dim> const>       mp_values;
    std::shared_ptr<BoundaryValues<dim> const> boundary_values;

    std::unordered_map<std::string, dealii::Vector<double> > vectors;
    std::unordered_map<std::string, double> values;
};

//////////////////////// SUPERCAPACITOR POSTPROCESSOR PARAMETERS ////////////////////////////
template <int dim>
class SuperCapacitorPostprocessorParameters : public PostprocessorParameters<dim> { 
public:
    SuperCapacitorPostprocessorParameters(std::shared_ptr<boost::property_tree::ptree const> d);
};

//////////////////////// SUPERCAPACITOR POSTPROCESSOR /////////////////////
template <int dim>
class SuperCapacitorPostprocessor : public Postprocessor<dim> {
public:
    SuperCapacitorPostprocessor(std::shared_ptr<PostprocessorParameters<dim> const> parameters);
    void reset(std::shared_ptr<PostprocessorParameters<dim> const> parameters);
private:
    bool debug_material_ids;
    bool debug_boundary_ids;
    std::vector<std::string> debug_material_properties;

};

} // end namespace cap

#endif // CAP_POSTPROCESSOR_H
