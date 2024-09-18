module CameraCalibrationMeta

using Rotations, CoordinateTransformations, StaticArrays, LinearAlgebra

export Calibration, lens_distortion, img2obj

const AM = AffineMap{Diagonal{Float64, SVector{2, Float64}}, SVector{2, Float64}}

const AMext = AffineMap{RotationVec{Float64}, SVector{3, Float64}}

const LM = LinearMap{SDiagonal{3, Float64}}

struct Calibration
    intrinsic::AM
    extrinsics::Vector{AMext}
    scale::LM
    k::Float64
    files::Vector{String}
    real2image
    image2real
end


"""
    lens_distortion
Lens distortion for one radial coefficient.
"""
function lens_distortion(v, k)
    k == 0 && return v
    r² = LinearAlgebra.norm_sqr(v)
    radial = 1 + k*r²
    radial*v
end

"""
    inv_lens_distortion
Analytical inverse lens distortion for one radial coefficient.
"""
function inv_lens_distortion(v2, k)
    k == 0 && return v2
    c = k*LinearAlgebra.norm_sqr(v2)
    rs = roots(Polynomial([c, 0, 1, -1]))
    rrs = filter(x -> abs(imag(x)) < 1e-10  , rs)
    radial = maximum(real, rrs)
    v2 / radial
end

# this is the inverse prespective map
depth(rc1, t, l) = -t/(l⋅rc1)
function get_inv_prespective_map(inv_extrinsic)
    function (rc)
        rc1 = push(rc, 1)
        t = inv_extrinsic.translation[3]
        l = inv_extrinsic.linear[end, :]
        d = depth(rc1, t, l)
        return d .* rc1
    end
end

function img2obj(intrinsic, extrinsics, scale, k)
    inv_extrinsics = inv.(extrinsics)
    inv_perspective_maps = get_inv_prespective_map.(inv_extrinsics)
    inv_distort(rc) = inv_lens_distortion(rc, k)
    return inv(scale), inv_extrinsics, inv_perspective_maps, inv_distort, inv(intrinsic)
end

end # module CameraCalibrationMeta
