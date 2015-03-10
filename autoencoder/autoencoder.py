from struct import *

def getTrainingLabels():
	f					= open("train-labels-idx1-ubyte", "rb")
	magic_number		= f.read(4)
	number_of_images	= f.read(4)
	noi					= unpack(">i",number_of_images)
	mn					= unpack(">i",magic_number)
	rest = unpack(">"+("B"*noi[0]),f.read(noi[0]))
	f.close()
	return noi[0],np.array(rest)	

def getTrainingImages():
	f					= open("train-images-idx3-ubyte","rb")
	mn					= unpack(">i",f.read(4))[0]
	noi					= unpack(">i",f.read(4))[0]
	rows				= unpack(">i",f.read(4))[0]
	cols				= unpack(">i",f.read(4))[0]
	
	images = []
	for i in range(noi):
			image = np.zeros((cols,rows))
			for r in range(rows):
				for c in range(cols):
					image[r,c]=unpack("B",f.read(1))[0]
			images.append(image)

	return noi,rows,cols,images
	
def getTestLabels():
	f					= open("t10k-labels-idx1-ubyte", "rb")
	magic_number		= f.read(4)
	number_of_images	= f.read(4)
	noi					= unpack(">i",number_of_images)
	mn					= unpack(">i",magic_number)
	rest = unpack(">"+("B"*noi[0]),f.read(noi[0]))
	f.close()
	return noi[0],np.array(rest)	

def getTestImages():
	f					= open("t10k-images-idx3-ubyte","rb")
	mn					= unpack(">i",f.read(4))[0]
	noi					= unpack(">i",f.read(4))[0]
	rows				= unpack(">i",f.read(4))[0]
	cols				= unpack(">i",f.read(4))[0]
	
	images = []
	for i in range(noi):
			image = np.zeros((cols,rows))
			for r in range(rows):
				for c in range(cols):
					image[r,c]=unpack("B",f.read(1))[0]
			images.append(image)

	return noi,rows,cols,images


##### define neurons #######

def logistic(z):
	return 1/(1+np.exp(-z))

def logistic_prime(a):
	return a*(1-a)

##### computer brain output ########

def computeActivation(inputs,weights,bias):
	return logistic(np.sum(inputs*weights)+bias)

def classify(image,weights_1,bias_1,weights_2,bias_2):
	num_inputs = len(image.flat)
	num_hiddens = len(bias_1.flat)
	num_outputs= len(bias_2.flat)

	activations_1 = np.zeros(num_hiddens)
	for i_hidden in range(num_hiddens):
		activations_1[i_hidden]=computeActivation(image.flatten(),weights_1[:,i_hidden].flatten(),bias_1[i_hidden])

	output = np.zeros(num_outputs)
	for i_output in range(num_outputs):
		output[i_output]=computeActivation(activations_1.flatten(),weights_2[:,i_output].flatten(),bias_2[i_output])

	return output,activations_1

##### train brain #######

def trainBrain(images,labels,diameter=28,hidden_diameter=24,iterations=-1,rate=0.01):
	weights_1 	= np.random.rand(diameter**2,hidden_diameter**2)/1000.
	bias_1 		= np.zeros(hidden_diameter**2)
	weights_2 	= np.random.rand(hidden_diameter**2,10)/1000.
	bias_2 		= np.zeros(10)

	if iterations == -1:
		iterations = len(images)

	label_list = []
	for i in range(10):
		label_list.append(np.zeros(10))
		label_list[i][i]=1

	sum_total = 0
	sum_correct = 0
	for i in range(iterations):
		output,activations_1 = classify(images[i]/256.,weights_1,bias_1,weights_2,bias_2)

		predicted = np.where(output==output.max())[0][0]
		sum_total+=1
		if predicted==labels[i]:
			sum_correct+=1
		print labels[i],predicted,output,sum_correct/float(sum_total)

		delta_weights_1 	= np.zeros((diameter**2,hidden_diameter**2))
		delta_bias_1 		= np.zeros(hidden_diameter**2)
		delta_weights_2 	= np.zeros((hidden_diameter**2,10))
		delta_bias_2 		= np.zeros(10)
		for j in range(10):
			cost_derivative			= output[j]-label_list[labels[i]][j]
			
			delta_bias_2[j]			+= cost_derivative * logistic_prime(output[j])
			delta_weights_2[:,j]	+= cost_derivative * logistic_prime(output[j]) * activations_1

			for k in range(hidden_diameter**2):	
				delta_bias_1[k]					+= cost_derivative * logistic_prime(output[j]) * weights_2[k,j] * logistic_prime(activations_1[k])
				delta_weights_1[:,k]			+= cost_derivative * logistic_prime(output[j]) * weights_2[k,j] * logistic_prime(activations_1[k]) * images[i].flatten()

		weights_1 	= weights_1 	- rate * delta_weights_1
		bias_1 		= bias_1 		- rate * delta_bias_1
		weights_2 	= weights_2 	- rate * delta_weights_2
		bias_2 		= bias_2		- rate * delta_bias_2

	return weights_1,bias_1,weights_2,bias_2


def applyModel(images,labels,weights_1,bias_1,weights_2,bias_2):
	iterations = len(images)

	label_list = []
	for i in range(10):
		label_list.append(np.zeros(10))
		label_list[i][i]=1

	sum_total = 0
	sum_correct = 0
	for i in range(iterations):
		output,activations_1 = classify(images[i]/256.,weights_1,bias_1,weights_2,bias_2)

		predicted = np.where(output==output.max())[0][0]
		sum_total+=1
		if predicted==labels[i]:
			sum_correct+=1
		print labels[i],predicted,output,sum_correct/float(sum_total)


###################### visualization routines ###########################################

def viewWeights(weights,inw,inh,outw,outh):
	pic = np.zeros((inh*outh,inw*outw))
	for j in range(outh):
		for i in range(outw):
			pic[j*inh:(j+1)*inh,i*inw:(i+1)*inw]=weights[:,j*outw+i].reshape(inh,inw)

	plt.imshow(pic,interpolation="nearest"); plt.colorbar(); plt.show()
