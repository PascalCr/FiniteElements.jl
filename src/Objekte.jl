export Platte

struct Platte{T<:AbstractFloat, Ti<:Integer}
    xL::T
    yL::T
    nx::Ti
    ny::Ti
    rho::Matrix{T}
    lambda::Matrix{T}
    cp::Matrix{T}
    temp::Matrix{T}
end


function Platte(xL::T, yL::T, nx::Ti, ny::Ti, dens, heat, cap, tv::Union{Number,Function}) where {T,Ti}
    if nx <= 0
        throw(ArgumentError("nx muss > 0 sein, ist aber $nx"))
    end

        function setgrid(f::Number)
            return [f for ix in 1:nx, iy in 1:ny]
        end

        function setgrid(f::Function)
            return [f((ix - 0.5) * xL / nx, (iy - 0.5) * xL / ny) for ix in 1:nx, iy in 1:ny]
        end

    rho = setgrid(dens)
    temp = setgrid(tv)
    lambda= setgrid(heat)
    cp = setgrid(cap)

    return Platte{T, Ti}(xL, yL, nx, ny, rho, lambda, cp, temp)
end

function tempiter(p::Platte, deltat::T) where {T}

    dtemp = [heattransfer(ix, iy, p) for ix in 1:p.nx, iy in 1:p.ny]
    tempn = dtemp * deltat + p.temp

    return Platte(p.xL, p.yL, p.nx, p.ny, p.rho, p.lambda, p.cp, tempn)
end

function heattransfer(ix, iy, p::Platte)

    produkt = [p.lambda[ix, iy] * nabla(p.temp, p)[ix, iy] for ix in 1:p.nx, iy in 1:p.ny]

    tempdiff = nabla(produkt, p)[ix, iy] / (p.cp[ix, iy] * p.rho[ix, iy])

end

function nabla(vm::Matrix{<:Vector}, p::Platte)
    vmx = [vm[ix, iy][1] for ix in 1:p.nx, iy in 1:p.ny]
    vmy = [vm[ix, iy][2] for ix in 1:p.nx, iy in 1:p.ny]
    [nablax(vmx, ix, iy, p) + nablay(vmy, ix, iy, p) for ix in 1:p.nx, iy in 1:p.ny]
end

function nabla(m::Matrix{<:Number}, p::Platte)
    [[nablax(m, ix, iy, p), nablay(m, ix, iy, p)] for ix in 1:p.nx, iy in 1:p.ny]
end

function nablax(m::Matrix, ix, iy, p::Platte)
    dx = p.xL / p.nx
    ix == 1 ? (m[ix + 1, iy] - m[ix, iy]) / dx :
    ix == p.nx ? (m[ix, iy] - m[ix - 1, iy]) / dx :
    (m[ix + 1, iy] - m[ix - 1, iy]) / (2 * dx)
end

function nablay(m::Matrix, ix, iy, p::Platte)
    dy = p.yL / p.ny
    iy == 1 ? (m[ix, iy + 1] - m[ix, iy]) / dy :
    iy == p.ny ? (m[ix, iy] - m[ix, iy - 1]) / dy :
    (m[ix, iy + 1] - m[ix, iy - 1]) / (2 * dy)
end
