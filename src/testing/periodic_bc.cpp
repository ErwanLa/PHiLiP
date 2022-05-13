#include <deal.II/base/tensor.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>


template <int dim>
void print_mesh_info(const Triangulation<dim> &triangulation,
                     const std::string &       filename)
{
  std::cout << "Mesh info:" << std::endl
            << " dimension: " << dim << std::endl
            << " no. of cells: " << triangulation.n_active_cells() << std::endl;
  {
    std::map<types::boundary_id, unsigned int> boundary_count;
    for (const auto &face : triangulation.active_face_iterators())
      if (face->at_boundary())
        boundary_count[face->boundary_id()]++;
    std::cout << " boundary indicators: ";
    for (const std::pair<const types::boundary_id, unsigned int> &pair :
         boundary_count)
      {
        std::cout << pair.first << '(' << pair.second << " times) ";
      }
    std::cout << std::endl;
  }
  std::ofstream out(filename);
  GridOut       grid_out;
  grid_out.write_vtu(triangulation, out);
  std::cout << " written to " << filename << std::endl << std::endl;
}


namespace PHiLiP {
namespace Tests {
template <int dim, int nstate>
/*AdvectionPeriodic<dim, nstate>::AdvectionPeriodic(const PHiLiP::Parameters::AllParameters *const parameters_input)
:
TestsBase::TestsBase(parameters_input)
{}*/


template <int dim, int nstate>
int AdvectionPeriodic<dim, nstate>::run_test() const
{
#if PHILIP_DIM==1
    using Triangulation = dealii::Triangulation<PHILIP_DIM>;
#else
    using Triangulation = dealii::parallel::distributed::Triangulation<dim>;
#endif

    std::shared_ptr<Triangulation> grid = std::make_shared<Triangulation>(
#if PHILIP_DIM!=1
        MPI_COMM_WORLD,
#endif
        typename dealii::Triangulation<dim>::MeshSmoothing(
            dealii::Triangulation<dim>::smoothing_on_refinement |
            dealii::Triangulation<dim>::smoothing_on_coarsening));

 double left = 0.0;
 double right = 2.0;
 const bool colorize = true;
 int n_refinements = 4;
 unsigned int poly_degree = 2;
 dealii::GridGenerator::hyper_cube(*grid, left, right, colorize);

 std::vector<dealii::GridTools::PeriodicFacePair<typename Triangulation::cell_iterator> > matched_pairs;
  dealii::GridTools::collect_periodic_faces(*grid,0,1,0,matched_pairs);
  dealii::GridTools::collect_periodic_faces(*grid,2,3,1,matched_pairs);
  //dealii::GridTools::collect_periodic_faces(*grid,4,5,2,matched_pairs);
  grid->add_periodicity(matched_pairs);

 grid->refine_global(n_refinements);

  dealii::GridTools::transform(
    [](const Point<dim> &in) {
      return Point<dim>(in[0]+ 0.5*std::sin(numbers::PI * in[1]), in[1] + 0.5*std::sin(numbers::PI * in[0] ));
    },
    grid);
  print_mesh_info(grid, "grid_periodic.vtu");


/* std::shared_ptr < PHiLiP::DGBase<dim, double> > dg = PHiLiP::DGFactory<dim,double>::create_discontinuous_galerkin(all_parameters, poly_degree, grid);
 dg->allocate_system ();*/

 return 0; //need to change
}

} //Tests namespace
} //PHiLiP namespace