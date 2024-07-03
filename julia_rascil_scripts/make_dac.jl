using DeconvMultiStep
using FITSIO

root = "."

lambda = parse(Float64, ARGS[1])
psf_high = read(FITS(joinpath(root, ARGS[2]))[1])
dirty_high = read(FITS(joinpath(root, ARGS[3]))[1])
psf_low = read(FITS(joinpath(root, ARGS[4]))[1])
dirty_low = read(FITS(joinpath(root, ARGS[5]))[1])
wavelet_dict = parse(Int, ARGS[6])
num_fista_iter = parse(Int, ARGS[7])
recon_noise = parse(Float64, ARGS[8])
vis_noise = parse(Float64, ARGS[9])

# make LP and HP filters
ℓ = parse(Float64, ARGS[10])
δ = parse(Float64, ARGS[11])
n_pix, _ = size(psf_high)

output_filename = joinpath(root, ARGS[12])

len = sqrt(recon_noise * recon_noise + vis_noise * vis_noise)

#recon_noise /= len
#vis_noise /= len

G = make_filters(ℓ, δ, n_pix, σ² = vis_noise, η² = recon_noise)

# reconstruction
if wavelet_dict == 0
    i_multistep = fista(psf_high, dirty_high, lambda, num_fista_iter; G=G, ip = dirty_low, Hl = psf_low)
else
    i_multistep, m = fista_iuwt(psf_high, dirty_high, lambda, num_fista_iter; G=G, ip = dirty_low, Hl = psf_low)
end

f = FITS(output_filename, "w")
write(f, i_multistep)
close(f)