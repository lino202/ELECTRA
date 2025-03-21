/*
 * ELECTRA. Electrophysiology Simulation Software.
 * Copyright (C) 2019
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/**
   \file cell_factory.hpp
   \brief CellFactory class header file.
   \author Konstantinos A. Mountris
   \date 23/10/2019
*/

#ifndef ELECTRA_ENGINE_PHYSICS_REACTION_DIFFUSION_FACTORY_HPP_
#define ELECTRA_ENGINE_PHYSICS_REACTION_DIFFUSION_FACTORY_HPP_

#include "ELECTRA/engine/physics/reaction_diffusion.hpp"
#include "ELECTRA/engine/physics/monodomain.hpp"
#include "ELECTRA/engine/physics/bidomain.hpp"
#include "ELECTRA/engine/utilities/logger.hpp"

#include <memory>
#include <string>


namespace ELECTRA {

/**
 *  \addtogroup Physics \{ */

/**
 * \class ReactionDiffusionFactory
 * \brief Class implementing a factory for reaction diffusion model generation.
 */
template <short DIM, short CELL_NODES>
class ReactionDiffusionFactory
{
public:

    /**
     * \brief Creates a reaction diffusion model according to the provided type.
     * \param [in] rd_model_type The type of the desired reaction diffusion model.
     * \return [std::shared_ptr<ReactionDiffusion>] Shared pointer to the created reaction diffusion model.
     */
    inline static std::shared_ptr<ReactionDiffusion<DIM,CELL_NODES>> Create(ReactionDiffusionType rd_model_type) {

        std::shared_ptr<ReactionDiffusion<DIM,CELL_NODES>> rd_ptr;

        switch (rd_model_type) {
        case ReactionDiffusionType::bidomain :
            //TODO Bidomain is deactivated as vout_ is deprecated now
            // for reactivating it, save the solution with the ensight exporter under 
            // Bidomain::Compute (similarly to how it is done for Monodomain::Compute) and test ulteriorly as this is a new feature
            throw std::invalid_argument(ELECTRA::Logger::Error("Bidomain is not supported for now"));
            // rd_ptr = std::make_shared<Bidomain<DIM,CELL_NODES>>();
            // break;
        case ReactionDiffusionType::monodomain :
            rd_ptr = std::make_shared<Monodomain<DIM,CELL_NODES>>();
            break;
        default:
            std::string error_msg = Logger::Error("Could not create reaction diffusion model. Not supported type.");
            throw std::invalid_argument(error_msg);
            break;
        }
        return rd_ptr;
    }
};


/** \} End of Doxygen Groups */

} // End of namespace ELECTRA

#endif // ELECTRA_ENGINE_PHYSICS_REACTION_DIFFUSION_FACTORY_HPP_