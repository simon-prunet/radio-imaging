using DeconvMultiStep
using FITSIO

root = "."

lambda = parse(Float64, ARGS[1])
psf = read(FITS(joinpath(root, ARGS[2]))[1])
dirty = read(FITS(joinpath(root, ARGS[3]))[1])
i_lowres = read(FITS(joinpath(root, ARGS[4]))[1])
wavelet_dict = parse(Int, ARGS[5])
num_fista_iter = parse(Int, ARGS[6])
recon_noise = parse(Float64, ARGS[7])
vis_noise = parse(Float64, ARGS[8])

# make LP and HP filters
ℓ = parse(Float64, ARGS[9])
δ = parse(Float64, ARGS[10])
n_pix, _ = size(psf)

output_filename = joinpath(root, ARGS[11])

len = sqrt(recon_noise * recon_noise + vis_noise * vis_noise)

recon_noise /= len
vis_noise /= len

G = make_filters(ℓ, δ, n_pix, σ² = vis_noise, η² = recon_noise)

# reconstruction
if wavelet_dict == 0
    i_multistep = fista(psf, dirty, lambda, num_fista_iter; G=G, ip = i_lowres)
else
    i_multistep, m = fista_iuwt(psf, dirty, lambda, num_fista_iter; G=G, ip = i_lowres)
end

f = FITS(output_filename, "w")
write(f, i_multistep)
close(f)