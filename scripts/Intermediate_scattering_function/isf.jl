
#Python version (isf_standalone_sinc.py) is faster than JL
start_time = time_ns()
using MDToolbox, LinearAlgebra, ArgParse, Printf, Base.Threads, Statistics, DelimitedFiles

t=mdload("asyn.psf")
t=mdload("asyn.dcd", top = t)

n_frames = t.nframe
println(n_frames)
const n_atoms = t["protein and atomname CA"].natom
println(n_atoms)

positions = t.xyz;

test_positions = positions[begin:10,:]
# k_vector in 3D
const k_vec = 1/7.2

function calculate_ISF(positions, max_t)
    lagtimes = 1:max_t
    timeseries = zeros(Float32, 2, max_t)
    timeseries[1,:] = lagtimes
    @threads for lag in lagtimes
        dr_vec = @view(positions[lag+1:end,:]) - @view(positions[begin:end-lag,:])
        dr_mul_k = dr_vec*k_vec
        # 
        mean_over_particles=zeros(Float32, size(dr_mul_k)[1])
        for i in 1:size(dr_mul_k)[1]
            temp_dr = reshape(@view(dr_mul_k[i,:]),(3,n_atoms))
            norm_array = @inbounds [norm(temp_dr[:,j]) for j in 1:size(temp_dr)[2]]
            # sinc_array = [sinc(x/pi) for x in norm_array]
            sinc_array = @fastmath @inbounds @. sinc(norm_array/pi)
            mean_over_particles[i] = @fastmath mean(sinc_array)
        end
        timeseries[2, lag] = @fastmath mean(mean_over_particles)
    end
    return timeseries
end

calculate_ISF(test_positions, 2)

@time res = calculate_ISF(positions, 500)

res_transpose = transpose(res)
writedlm( "ISF_JULIA.csv",  res_transpose, '\t')
# println(res)

end_time = time_ns()
dt = (end_time - start_time) / 10^9
println("Execution time: ", dt, "(s)")
