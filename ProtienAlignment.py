"""
Aligning 2 protien sequences using dynamic Programming (Needleman Wunsch)
using identity scoring scheme(no gap penalty)
Output: Alignement of protein sequences and 
Storing the dotplot, sum matrix in csv files(',' separated)


"""
import pandas as pd
import copy
import sys,getopt

def main(argv):

	inputfile = 'protein.fa'
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	except getopt.GetoptError:
		print ('Run the file in this format "python3','test.py -i <inputfile> "')
		sys.exit(1)
	for opt, arg in opts:
		if opt in ("-i", "--ifile"):
			inputfile = arg

			
	#Step0: extracting data of 2 protien sequences from a fasta file

	seq=[]
	w=open(inputfile)
	for line in w:
		if(line[:1]!='>' and line!='\n'):
			seq.append(line[:-1]);

	print("PROTIEN SEQUENCES BEFORE ALLIGNMENT\n")	
	print(seq[0])
	print(seq[1],"\n")


	#Step1: Making a dotplot of the 2 protien sequences

	w=len(seq[0])
	h=len(seq[1])
	dotplot = [[0 for x in range(w)] for y in range(h)]  

	for j in range(len(seq[0])):
		for k in range(len(seq[1])):
			if(seq[0][j]==seq[1][k]):
				dotplot[j][k]=1

	# making a visual dataframe and converting to a csv file(',' separated) 

	colnames=[]
	rownames=[]

	for c in seq[0][:]:
		colnames.append(c)
	for r in seq[1][:]:
		rownames.append(r)

	df=pd.DataFrame(dotplot,index=[colnames],columns=[x for x in rownames])
	df=(df.transpose())
	print(df,'\n')
	df.to_csv('dotplot_output.csv',sep=',',encoding='utf-8')


	#Step2 : Making a sum matrix and storing it in a csv file(',' separated)


	sum_matrix=copy.deepcopy(dotplot)

	for j in range(len(seq[0])-2,-1,-1):
		for k in range(len(seq[1])-2,-1,-1):
			if(j+2<h):
				max_colval=0
				for p in range(j+2,len(seq[1])-1):
					if(sum_matrix[p][k+1]>max_colval):
						max_colval=sum_matrix[p][k+1]
			if( k+2<w and j+2<h):
				sum_matrix[j][k]+=max(sum_matrix[j+1][k+1],max(sum_matrix[j+1][k+2:]),max_colval)
				
			elif (k+2<w):
				sum_matrix[j][k]+=max(sum_matrix[j+1][k+1],max(sum_matrix[j+1][k+2:]))
			
			elif (j+2<h):
				sum_matrix[j][k]+=max(sum_matrix[j+1][k+1],max_colval)

			else:
				sum_matrix[j][k]+= (sum_matrix[j+1][k+1])

	df_sum=pd.DataFrame(sum_matrix,index=[colnames],columns=[x for x in rownames])
	df_sum=(df_sum.transpose())
	print(df_sum)
	df_sum.to_csv('sum_matrix_output.csv')

	#Step3: backtracing and finding the required allignment by introducing gaps

	diag=0
	first_seq_allign=''
	sec_seq_allign=''
	max_index=0
	prevmax_index=-1
	for j in range(len(sum_matrix[0])):
		max_diag=0
		
		for k in range(diag,len(sum_matrix[0])):
			if sum_matrix[j][k]>max_diag:
				max_diag=sum_matrix[j][k]
				max_index=k


		if(prevmax_index-max_index==0):
			#introduce gaps in 2nd sequence
			sec_seq_allign+='_'
			first_seq_allign+=seq[0][j]
			diag-=1

		elif(j<max_index and dotplot[j][max_index]==1):
			#introduce gaps in 1st seq
			ngaps=(max_index - prevmax_index)-1

			for g in range(ngaps):
				first_seq_allign+='_'
				sec_seq_allign+=seq[1][prevmax_index+g+1]
			diag+=ngaps
			first_seq_allign+=seq[0][j]
			sec_seq_allign+=seq[1][max_index]

		elif(j<max_index):
			first_seq_allign+=seq[0][j]
			sec_seq_allign+=seq[1][max_index]
		
		elif(j>max_index and dotplot[j][max_index]==1):
			first_seq_allign+=seq[0][j]
			sec_seq_allign+=seq[1][max_index]

		elif(sum_matrix[j][max_index]==0):
			sec_seq_allign+='_'
			first_seq_allign+=seq[0][j]
		else:
			first_seq_allign+=seq[0][j]
			sec_seq_allign+=seq[1][max_index]

		prevmax_index=max_index
		diag+=1


	#The final best possible allignments of 2 sequences are:
	print("\nPROTIEN SEQUENCES AFTER ALLIGNMENT\n")
	print(first_seq_allign)
	print(sec_seq_allign)


if __name__ == "__main__":
	main(sys.argv[1:])