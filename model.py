import keras
from compoundPrep import MoleculePrepper
import numpy as np
from pharmaDb import engine, diseases, maxVecNormTbl

class PredictiveModel:

	#As it relates to compounds, only SMILES will be input here.

	#TODO: mol2vec functionality must be provided.
	#TODO: Design db and take into consideration the disease representations.

	def __init__(self, encoderFileName=" ", predictiveModelName=" "):
		#Load encoder model.
		self.encoder = keras.models.load_model(encoderFileName)
		#Load predictive model.
		self.predictiveModel = keras.models.load_model(predictiveModelName)
		

		self.resultsDict = dict()

		self.molVector = None
		self.SMILES = None
		self.diseaseVector = None
		self.absorption = None
		self.distribution = None
		self.metabolism = None
		self.excretion = None
		self.toxicity = None
		self.xlogP = None
		self.pKa = None
		self.solubility = None
		self.size = None
		self.stability = None
		self.molecularWeight = None
		self.mass = None
		self.efficacy = None

		#A dictionary to store the disease name and corresponding vector representation.
		self.diseaseDict = dict()
		#TODO: Fill out this dictionary.
		#NOTE: The below is temporary. We must 

		#Note: we must be pulling this data from the database.
		conn = engine.connect()
		diseasesList = conn.execute(diseases.select())
		for row in diseasesList:
			self.diseaseDict[row.id] = self.__vecStrToList(row.diseaseRepresentation)


		#NOTE: This list will specify what features the model has been trained on and the order in which 
			#the features must be specified in the input list.
		#TODO: Set the order which the features must be fed into the model.
		'''
		Note that we will only be including the features that are used for the current model.
		'''
		#self.featuresOrdering = ['molVector', 'diseaseVector', 'absorption', 'distribution', 'metabolism', 'excretion','toxicity', 'xlogP','pKa','solubility','size','stability','molecularWeight', 'mass']
		self.featuresOrdering = ['molVector', 'diseaseVector','xlogP', 'molecularWeight', 'mass']
		
		set(['SMILES', 'diseaseVector', 'absorption', 'distribution', 'metabolism', 'excretion','toxicity', 'xlogP','pKa','solubility','size','stability','molecularWeight', 'mass'])



		self.modelInputList = []
		self.compressedModelInputList = []

		
		#To take care of preparing the input molecule.
		self.molPrepper = MoleculePrepper()


		#Dictionary to store the 
		self.funcDict = dict()

		self.funcDict['SMILES'] = self.molPrepper.getSmiles
		self.funcDict['absorption'] = self.molPrepper.getAbsorption
		self.funcDict['distribution'] = self.molPrepper.getDistribution
		self.funcDict['metabolism'] = self.molPrepper.getMetabolism
		self.funcDict['excretion'] = self.molPrepper.getExcretion
		self.funcDict['toxicity'] = self.molPrepper.getToxicity
		self.funcDict['xlogP'] = self.molPrepper.getXLogP
		self.funcDict['pKa'] = self.molPrepper.getPKA
		self.funcDict['solubility'] = self.molPrepper.getSolubility
		self.funcDict['size'] = self.molPrepper.getSize
		self.funcDict['molecularWeight'] = self.molPrepper.getWeight
		self.funcDict['stability'] = self.molPrepper.getStability
		self.funcDict['mass'] = self.molPrepper.getMass
		
		conn = engine.connect()
		tbl = conn.execute(maxVecNormTbl.select())
		self.maxVecNorm = None
		for row in tbl:
			self.maxVecNorm = row.maxVecNorm
			break

	def reset(self):
		self.molPrepper.reset()
	'''
	Used for setting the features that will be used in the model.
	Not sure if we will use it yet.
	'''
	def setModelFeatures(self, featuresSet):
		pass





	def __vecStrToList(self, vecStr):
		tempList = []

		if vecStr == None:
			return None

		vecList = vecStr.split(',')
		#print(len(vecList))
		for numStr in vecList:
			if numStr == '' or numStr == ' ':
				continue
			tempList.append(float(numStr))

		return tempList



	def displayFeaturesModelTrainedOn(self):
		pass


	def displayKeys(self):
		pass

	#Returns a list.
	def __smilesToVec(self, smiles):
		pass

	#Check to determine if the input is compatible.
	def __checkCompatability(self, inputFeatures):
		pass


	def __compressData(self, inputArray):
		#print(inputArray.shape)
		return self.encoder.predict(np.array([inputArray]))[0]


	#Input vector must be input as an np.array object [...,...]

	'''
	We will allow the user to enter the molecule name or structure. But for now we will only have them enter
	the molecule name. We need to expand the functionality to include being able to make predictions
	given the molecular structure.
	'''
	def predict(self, molName=None, molStructure=None):
		#Retrieve the properties given the molecule name.
		#We will assume for now the user is only going to be entering the molecule name.

		#Note that a molecule may have several/aliases, so we need to take care of this.
		#Get smiles.
		self.molPrepper.setMoleculeFromName(molName)

		molPropertiesDict = dict()

		#Get molecule vector.
		v = self.molPrepper.molToVec()
		#print(v.shape)
		
		print(type(v) == np.ndarray) 
		if type(v) != np.ndarray:
			self.resultsDict = None
			return None

		print(type(v))

		molPropertiesDict['molVector'] = np.array(self.__compressData(v))
		#Now we need to normalize the vector.
		tempVec = []
		for elmt in molPropertiesDict['molVector']:
			tempVec.append(elmt/self.maxVecNorm)

		molPropertiesDict['molVector'] = np.array(tempVec)

		#Compute properties.
		featuresOrderingTemp = self.featuresOrdering.copy()	#Copy the elements.
		featuresOrderingTemp.remove('molVector')
		featuresOrderingTemp.remove('diseaseVector')

		for feature in featuresOrderingTemp:
			featureVal = self.funcDict[feature]() #Invoke function to compute property.
			print(feature)
			if featureVal == None:
				self.resultsDict = None
				return None	#We need all of the features requested or else we cannot make predictions.
			molPropertiesDict[feature] = featureVal



		#TODO Describe what this is for.
		self.resultsDict = dict()

		'''
		Remember that we have to feed the data into the model for each disease.
		This involves preparing the data inputs for each disease.
		'''
		for diseaseKey in self.diseaseDict:
			diseaseVector = self.diseaseDict[diseaseKey]	#Disease represented as a vector (list).

			dataInput = []

			#Remember that we must flatten the data; that is, a will not be storing a list within a list.
			#Also remember that the input form is first molecule vector then disease vector,....
			
			for molVecElmt in molPropertiesDict['molVector']:
				dataInput.append(molVecElmt)

			for diseaseVecElmt in diseaseVector:
				dataInput.append(diseaseVecElmt)

			#Now iterate over the remaining feature set and add them to the input list in the correct order.
			#Also note that the rest of the inputs are scalars, so we do not have to flatten any lists.
			for feature in featuresOrderingTemp:
				dataInput.append(molPropertiesDict[feature])
				#print(feature)
				#print(molPropertiesDict[feature])
				#print("\n\n")


			dataInput = np.array([dataInput])

			modelOutput = self.predictiveModel.predict(dataInput)[0]

			#Store this in a dictionary.
			#Key = disease key, value = predicted effectiveness (in arange [0,1]).
			self.resultsDict[diseaseKey] = modelOutput

	#Returns the results dictionary. Users can iterate over this and get 
	def getResultsDict(self):
		return self.resultsDict









