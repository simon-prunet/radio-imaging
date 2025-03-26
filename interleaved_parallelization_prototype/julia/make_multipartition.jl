using DeconvMultiStep
using FITSIO
using FFTW

function imfilter(img::Matrix{T}, psf::Matrix{T}) where {T<:Float64}

    @assert size(img) == size(psf) "Image and psf must have the same size"
    real(ifft(fft(img).*fft(ifftshift(psf))))
end

#for now read from command line, if this poses a problem, change to reading from JSON file

root = "."

psf = read(FITS(joinpath(root, ARGS[1]))[1])
dirty = read(FITS(joinpath(root, ARGS[2]))[1])
lambda = parse(Float64, ARGS[3])
n_fista_iter = parse(Int64, ARGS[4])
num_partitions = parse(Int64, ARGS[5])
curr_partition = parse(Int64, ARGS[6])
curr_maj_cycle = parse(Int64, ARGS[7])
delta = parse(Int64, ARGS[8])
output_filename = ARGS[9]

deconvolved = Nothing

if curr_maj_cycle > 0
	curr_arg_idx = 9
	curr_constr_idx = 1
	constraint_images = Matrix{Float64}[]
	for i in 1:num_partitions
		push!(constraint_images, read(FITS(joinpath(root, ARGS[curr_arg_idx + curr_constr_idx]))[1]))
		global curr_constr_idx += 1
	end

	curr_arg_idx += num_partitions
	ells = Int64[]
	for i in 1:num_partitions - 1
		push!(ells, parse(Int64, ARGS[curr_arg_idx + i]))
	end

	curr_arg_idx += num_partitions - 1
	sigma2s = Float64[]
	for i in 1:num_partitions
		push!(sigma2s, parse(Float64, ARGS[curr_arg_idx + i]))
	end

	n_pix, _ = size(psf)

	G = make_filters_multipartition(ells, delta, sigma2s, n_pix) 

	deconvolved = fista_multipartition(psf, dirty, lambda, n_fista_iter, G, constraint_images, curr_partition)
else
	deconvolved = fista(psf, dirty, lambda, n_fista_iter)
end

f = FITS(joinpath(root, output_filename), "w")
write(f, deconvolved)
close(f)