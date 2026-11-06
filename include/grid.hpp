#pragma once

#include <cstddef>

namespace ccp {

/// Uniform 1D grid, staggered in the usual finite-volume way.
///
///   cells   |  0  |  1  |  2  | ... | N-1 |     densities, potential, energy
///   faces   0     1     2     3         N-1  N   fluxes, electric field
///
/// The electrodes sit on faces 0 and N. Densities live at cell centres so that
/// the divergence of a face flux is a clean difference; the field lives on
/// faces so drift velocities are defined exactly where the fluxes need them.
class Grid {
public:
    Grid(std::size_t num_cells, double length)
        : num_cells_(num_cells), length_(length),
          dx_(length / static_cast<double>(num_cells)) {}

    std::size_t num_cells() const noexcept { return num_cells_; }
    std::size_t num_faces() const noexcept { return num_cells_ + 1; }
    double length() const noexcept { return length_; }
    double dx() const noexcept { return dx_; }

    /// Position of the centre of cell i.
    double cell_centre(std::size_t i) const noexcept {
        return (static_cast<double>(i) + 0.5) * dx_;
    }

    /// Position of face f (face 0 is the powered electrode, face N the ground).
    double face_position(std::size_t f) const noexcept {
        return static_cast<double>(f) * dx_;
    }

private:
    std::size_t num_cells_;
    double length_;
    double dx_;
};

}  // namespace ccp
