# Parameters for the solver
en_vocab=$1
fr_vocab=$2
dense_en_vectors=$3
dense_fr_vectors=$4
alignment_matrix=$5

# Hyperparameters
m_en=0.5
m_fr=0.5
m_xling=10
size=100
tol=1e-07		# Tolerance 
iters=10		# Num of iterations to run - 1000 is a good idea


matlab -nosplash -nojvm -nodisplay -r \
	"fasta_biling_solver(${m_en},${m_fr},${m_xling},${size},${tol},${iters},'${dense_fr_vectors}','${dense_en_vectors}','${alignment_matrix}')"

paste "${fr_vocab}" bisparse_fr.txt > final_bisparse_fr.txt
paste "${en_vocab}" bisparse_en.txt > final_bisparse_en.txt
rm bisparse*.txt
