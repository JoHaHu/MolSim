#include "XMLFileReader.h"
#include "SimulationInputSchema.hxx"
#include "container/boundary.h"
#include <iostream>
#include <memory>
#include <spdlog/spdlog.h>
#include <vector>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/PlatformUtils.hpp>

using namespace xercesc;

/**
 * @brief Parses the data from the given XML file using the XSD library and stores it in a Config object.
 *
 * @param xmlFilePath path of the XML file that contains the simulation parameter values
 *
 * creates a Config object to store the simulation parameters given by the XML input
 * for further passing and use inside the project as well as easier testing
 *
 */
auto XMLFileReader::parseXMLData(const std::string &xmlFilePath) -> std::shared_ptr<config::Config> {
    try {
        XMLPlatformUtils::Initialize();
    } catch (const XMLException &exception) {
        char *message = XMLString::transcode(exception.getMessage());
        SPDLOG_WARN("Error during XML file initialization!\n");
        SPDLOG_INFO(message);
        return nullptr;
    }

    try {
        // Parse the XML file and disable validation (according to worksheet)
        std::unique_ptr<scenario> scenario = scenario_(xmlFilePath, 0);
        auto &sc = *scenario;
        SPDLOG_INFO("Reading XML file.");

        // Log the header input values so the user can confirm their correctness
        SPDLOG_DEBUG("Passing the following header parameters to simulation: " + sc.header().base_name() +
                     " | output frequency: " + std::to_string(sc.header().output_frequency()) + " | t_end value: " +
                     std::to_string(sc.header().t_end()));

        config::Config config = config::Config();

        // Header parsing
        config.base_name = sc.header().base_name();
        config.end_time = sc.header().t_end();
        config.output_frequency = sc.header().output_frequency();
        config.output_file = sc.header().output_file();
        config.seed = sc.header().seed();
        config.delta_t = sc.header().delta_t();
        config.parallelized = sc.header().parallelized();
        config.vectorized = sc.header().vectorized();
        config.dimensions = sc.header().dimensions();

        /** Particle Loader choice **/

        // Linked Cells
        if (sc.container().linked_cells().present()) {
            config.particle_container_type = ParticleContainerType::LinkedCells;

            // Parsing array of domain size (volume)
            const auto &domain_size = sc.container().linked_cells()->domain_size();
            std::vector<double> domain_size_vector = {domain_size.x(), domain_size.y(), *domain_size.z()};

            config.domain_size = domain_size_vector;

            // Parsing array of domain size (volume)
            const auto &boundary_conditions = sc.container().linked_cells()->boundary_conditions();
            std::vector<BoundaryCondition> boundary_conditions_vec;
            for (const auto &arr: boundary_conditions.boundary_condition()) {
                BoundaryCondition bc = BoundaryCondition::none;
                if (arr.type() == "outflow") {
                    bc = BoundaryCondition::outflow;
                } else if (arr.type() == "reflecting") {
                    bc = BoundaryCondition::reflecting;
                } else if (arr.type() == "periodic") {
                    bc = BoundaryCondition::periodic;
                } else {
                    spdlog::warn("unknown boundary condition {}", arr.type());
                }
                boundary_conditions_vec.push_back(bc);
            }

            config.boundary_conditions = boundary_conditions_vec;
            config.cutoff_radius = sc.container().linked_cells()->cutoff_radius();

            // Vector
        } else if (sc.container().vector().present()) {
            config.particle_container_type = ParticleContainerType::Vector;
        }

        // Gravity Model
        //    if (sc.forces().gravity().present()) {
        //
        //      config.simulation_type = simulator::physics::ForceModel::Gravity;
        //
        //      std::vector<CelestialBody> temp_bodies;
        //      // Iterate through all celestial bodies
        //      const auto &bodies = sc.forces().gravity()->celestial_body();
        //      int body_count = 0;
        //      for (const auto &body : bodies) {
        //
        //        std::vector<double> coordinate{};
        //        std::vector<double> velocity{};
        //        double mass = body.mass();
        //
        //        for (const auto &arr : body.coordinate().value()) {
        //          coordinate.push_back(arr);
        //        }
        //
        //        for (const auto &arr : body.velocity().value()) {
        //          velocity.push_back(arr);
        //        }
        //
        //        // initializing arrays for the config object to be created
        //        std::array<double, 3> arr_coordinate{};
        //        std::array<double, 3> arr_velocity{};
        //
        //        // checking if lengths match before copying to static length arrays of size 3, otherwise give error
        //        if (coordinate.size() == arr_coordinate.size() && velocity.size() == arr_velocity.size()) {
        //          std::copy(coordinate.begin(), coordinate.end(), arr_coordinate.begin());
        //          std::copy(velocity.begin(), velocity.end(), arr_velocity.begin());
        //        } else {
        //          SPDLOG_WARN("There was an error while parsing the values for the cuboids.");
        //          return nullptr;
        //        }
        //
        //        auto temp_body = CelestialBody(arr_coordinate, arr_velocity, mass);
        //        temp_bodies.emplace_back(temp_body);
        //        body_count++;
        //      }
        //
        //      config.total_bodies = body_count;
        //      config.celestial_bodies = temp_bodies;
        //    }

        // Lennard-Jones Model
        if (sc.forces().lennard_jones().present()) {
            auto ljf = *sc.forces().lennard_jones();

            std::vector<double> sigma{};
            std::vector<double> epsilon{};
            std::vector<double> mass{};

            for (auto &type: ljf.particleTypes().particleType()) {
                sigma.emplace_back(type.sigma());
                epsilon.emplace_back(type.epsilon());
                mass.emplace_back(type.mass());
            }
            config.sigma = sigma;
            config.epsilon = epsilon;
            config.mass = mass;

            config.ljf_gravity = *ljf.gravity();

            // parsing settings of Lennard Jones force simulation model
            config.simulation_type = simulator::physics::ForceModel::LennardJones;

            if (ljf.particles().present()) {

                // Cuboids
                std::vector<Cuboid> temp_cuboids;
                // Iterate through all cuboid elements
                for (const auto &cuboid: ljf.particles()->cuboid()) {
                    std::vector<double> coordinate{cuboid.coordinate().x(), cuboid.coordinate().y(),
                                                   *cuboid.coordinate().z()};
                    std::vector<double> dimensions{cuboid.dimensions().x(), cuboid.dimensions().y(),
                                                   *cuboid.dimensions().z()};
                    std::vector<double> velocity{cuboid.velocity().x(), cuboid.velocity().y(), *cuboid.velocity().z()};


                    if (cuboid.fixed().present()) {
                        int flag = *cuboid.fixed();
                        if (flag == 1) {
                            auto temp_cuboid = Cuboid(coordinate, dimensions, velocity,
                                                      mass[cuboid.particleTypeId()], *cuboid.spacing(),
                                                      static_cast<int>(cuboid.particleTypeId()), 1);
                            temp_cuboids.emplace_back(temp_cuboid);
                        }
                    } else {
                        auto temp_cuboid = Cuboid(coordinate, dimensions, velocity,
                                                  mass[cuboid.particleTypeId()], *cuboid.spacing(),
                                                  static_cast<int>(cuboid.particleTypeId()), 0);
                        temp_cuboids.emplace_back(temp_cuboid);
                    }
                }
                config.cuboids = temp_cuboids;

                // Discs
                std::vector<Disc> temp_discs;

                // Iterate through all disc elements
                for (const auto &disc: ljf.particles()->disc()) {
                    std::vector<double> coordinate = {disc.coordinate().x(), disc.coordinate().y(),
                                                      *disc.coordinate().z()};
                    std::vector<double> velocity = {disc.velocity().x(), disc.velocity().y(), *disc.velocity().z()};
                    int radius = disc.radius();
                    auto temp_disc = Disc(coordinate, velocity, radius);
                    temp_discs.emplace_back(temp_disc);
                }
                config.discs = temp_discs;

                // Spheres

                std::vector<Sphere> temp_spheres;
                // Iterate through all disc elements
                for (const auto &sphere: ljf.particles()->sphere()) {

                    std::vector<double> coordinate = {sphere.coordinate().x(), sphere.coordinate().y(),
                                                      *sphere.coordinate().z()};
                    std::vector<double> velocity = {sphere.velocity().x(), sphere.velocity().y(),
                                                    *sphere.velocity().z()};
                    int radius = sphere.radius();

                    auto temp_sphere = Sphere(coordinate, velocity, radius, sphere.mesh_width());
                    temp_spheres.emplace_back(temp_sphere);
                }

                config.spheres = temp_spheres;

                std::vector<Torus> temp_tori;
                // Iterate through all disc elements
                for (const auto &torus: ljf.particles()->torus()) {

                    std::vector<double> coordinate = {torus.coordinate().x(), torus.coordinate().y(),
                                                      *torus.coordinate().z()};
                    std::vector<double> velocity = {torus.velocity().x(), torus.velocity().y(), *torus.velocity().z()};

                    double major_radius = torus.major_radius();
                    double minor_radius = torus.minor_radius();
                    auto temp_torus = Torus(coordinate, velocity, major_radius, minor_radius);
                    temp_tori.emplace_back(temp_torus);
                }
                config.tori = temp_tori;

                std::vector<DoubleHelix> temp_helices;
                // Iterate through all disc elements
                for (const auto &helix: ljf.particles()->doubleHelix()) {

                    std::vector<double> coordinate = {helix.coordinate().x(), helix.coordinate().y(),
                                                      *helix.coordinate().z()};
                    std::vector<double> velocity = {helix.velocity().x(), helix.velocity().y(), *helix.velocity().z()};
                    double radius = helix.radius();
                    double pitch = helix.pitch();
                    double height = helix.height();
                    auto temp_helice = DoubleHelix(coordinate, velocity, radius, pitch, height);
                    temp_helices.emplace_back(temp_helice);
                }
                config.double_helices = temp_helices;
            }
        }

        if (scenario->thermostat().present()) {
            auto thermostat = *scenario->thermostat();

            config.temp_init = thermostat.t_init();
            config.temp_target = thermostat.t_target().present() ? *thermostat.t_target() : std::optional<double>();
            config.max_temp_diff = thermostat.max_temp_diff().present() ? *thermostat.max_temp_diff()
                                                                        : std::optional<double>();
            config.thermo_step = thermostat.frequency();
            config.use_brownian_motion = thermostat.brownian_motion().present();
            config.brownian_motion = *thermostat.brownian_motion();
        } else {
            config.thermo_step = 0;
        }

        if (scenario->checkpoints().present()) {
            auto input_checkpoints = std::vector<std::string>();
            if (scenario->checkpoints()->path().present()) {
                config.output_checkpoint = *scenario->checkpoints()->path();
            }

            auto cps = *scenario->checkpoints();
            for (auto cp: cps.checkpoint()) {
                input_checkpoints.push_back(cp);
            }
            config.input_checkpoints = input_checkpoints;
        }

        return std::make_shared<config::Config>(config);

    } catch (const xml_schema::exception &exception) {
        SPDLOG_WARN("An error has occurred! Please look at the exception details here:\n");
        std::cerr << exception << '\n';
    }
    XMLPlatformUtils::Terminate();
    return nullptr;
}