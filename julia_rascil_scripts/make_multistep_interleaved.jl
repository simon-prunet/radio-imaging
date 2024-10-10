using DeconvMultiStep
using FITSIO

root = "."

lambda = parse(Float64, ARGS[1])
psf = read(FITS(joinpath(root, ARGS[2]))[1])
dirty = read(FITS(joinpath(root, ARGS[3]))[1])
i_constraint = read(FITS(joinpath(root, ARGS[4]))[1])
wavelet_dict = parse(Int, ARGS[5])
num_fista_iter = parse(Int, ARGS[6])
vl_variance = parse(Float64, ARGS[7])
vh_variance = parse(Float64, ARGS[8])

# make LP and HP filters
ℓ = parse(Float64, ARGS[9])
δ = parse(Float64, ARGS[10])
n_pix, _ = size(psf)

partition = parse(Int, ARGS[11])
maj_cycle = parse(Int, ARGS[12])

output_filename = joinpath(root, ARGS[13])



# reconstruction
if maj_cycle == 0
    if wavelet_dict == 0
        i_multistep = fista(psf, dirty, lambda, num_fista_iter)
    else
        i_multistep, m = fista_iuwt(psf, dirty, lambda, num_fista_iter)
    end
else
    len = sqrt(vl_variance * vl_variance + vh_variance * vh_variance)

    vl_variance /= len
    vh_variance /= len

    switch = partition == 0

    G = make_filters(ℓ, δ, n_pix, σ² = vh_variance, η² = vl_variance, switch = switch)

    if wavelet_dict == 0
        i_multistep = fista(psf, dirty, lambda, num_fista_iter; G=G, ip = i_constraint)
    else
        i_multistep, m = fista_iuwt(psf, dirty, lambda, num_fista_iter; G=G, ip = i_constraint)
    end
end

f = FITS(output_filename, "w")
write(f, i_multistep)
close(f)