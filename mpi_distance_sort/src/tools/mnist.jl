using LinearAlgebra, MLDatasets

"""
Write a matrix into a binary file (to read from C)
"""
function write_matrix_binary(A::AbstractMatrix{Float64}, filename::String)

	open(filename, "w") do file
    	write(file, size(A,1)) # write an Int64
    	write(file, size(A,2)) # write an Int64
		for i in eachindex(A)
			write( file, A[i] ) # write Float64 (column-major by default)
		end
    end
end

train_x, _ = MNIST.traindata()
test_x,  _ = MNIST.testdata()
X = Float64.( [reshape(train_x,28*28,:) reshape(test_x,28*28,:)] )
write_matrix_binary(X, "/home/vasilisk/Git/pds_assignment_2/mpi_distance_sort/datasets/mnist.bin")