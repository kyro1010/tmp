#Molecule preparation
#This will provide functionality to convert a molecule name to a smiles representation.

#We need to allow the user to enter a molecular structure and we need to be able to conver it to SMILES.

#We need to allow the user to enter the molecule name and we need to convert it to SMILES.

#Check out pubChemPy and rdkit.

import pubchempy as pcp
import numpy as np
import pandas as pd
import os
from tqdm import tqdm_notebook as tqdm

from rdkit import Chem
from mol2vec.features import mol2alt_sentence
import pybel

from mol2vec.features import mol2alt_sentence, mol2sentence, MolSentence, DfVec, sentences2vec

from gensim.models import Word2Vec

from pharmaDb import engine, drugsOriginal, featureMaps


class MoleculePrepper:

	def __init__(self):
		self.molecule = None
		self.smiles = None
		self.vector = None
		self.model = Word2Vec.load('model_300dim.pkl')

		#self.cannotUseInPredictiveModel = False

		#We need to query the database to get the maps we need.
		conn = engine.connect()
		fm = conn.execute(featureMaps.select()).fetchone()
		
		self.mapDict = dict()
		self.mapDict['mass'] = dict()
		self.mapDict['xlogP'] = dict()
		self.mapDict['molecular_weight'] = dict()

		if fm != None:
			#print("\n\n\nfm\n\n\n")
			self.mapDict['xlogP']['min'] = fm.xlogP_min
			self.mapDict['xlogP']['max'] = fm.xlogP_max
			self.mapDict['molecular_weight']['min'] = fm.molecularWeight_min
			self.mapDict['molecular_weight']['max'] = fm.molecularWeight_max
			self.mapDict['mass']['min'] = fm.mass_min
			self.mapDict['mass']['max'] = fm.mass_max

		'''
		We will no longer be making API calls to PubChem from here. All of our calls from here will be
		to our own database.
		Query database.
		Store results in a dictionary.
		'''

		self.drugsDict = dict()	#A dictionary that is used for storing the properties of various drugs in our system.
		drugsTbl = conn.execute(drugsOriginal.select())
		for row in drugsTbl:
			drugName = row.drugName
			self.drugsDict[drugName] = dict()
			self.drugsDict[drugName]['SMILES'] = row.SMILES
			self.drugsDict[drugName]['vector'] = row.vector
			self.drugsDict[drugName]['absorption'] = row.base_absorption
			self.drugsDict[drugName]['distribution'] = row.base_distribution
			self.drugsDict[drugName]['metabolism'] = row.base_metabolism
			self.drugsDict[drugName]['excretion'] = row.base_excretion
			self.drugsDict[drugName]['toxicity'] = row.base_toxicity
			self.drugsDict[drugName]['xlogP'] = row.base_xlogP
			self.drugsDict[drugName]['solubility'] = row.base_solubility
			self.drugsDict[drugName]['size'] = row.base_size
			self.drugsDict[drugName]['stability'] = row.base_stability
			self.drugsDict[drugName]['molecular_weight'] = row.base_molecularWeight
			self.drugsDict[drugName]['mass'] = row.base_mass



	def reset(self):
		self.molecule = None
		self.smiles = None
		self.vector = None



	def setMoleculeFromName(self, molName):
		#self.molecule = pcp.get_compounds(molName,"name")[0]
		self.molecule = molName #We will no longer be making API calls to pubchem from here. We will only be making calls to our own database.


	def setMoleculeFromSmiles(self, inputSmiles):
		#self.molecule = #Get the molecule name from the SMILES representation.
		self.smiles = inputSmiles

	#Method to set the molecule from the molecule structure.
	def setMoleculeFromStructure(self, molStructure):
		pass

	def getAbsorption(self):
		pass

	def getDistribution(self):
		pass

	def getMetabolism(self):
		pass

	def getExcretion(self):
		pass

	def getToxicity(self):
		pass

	def getXLogP(self):		
		'''xlogP = self.molecule.xlogp
		if xlogP == None:
			self.cannotUseInPredictiveModel = True
			return None'''
		#xlogP = float(self.molecule.xlogp)
		if self.molecule not in self.drugsDict:
			return None

		xlogP = self.drugsDict[self.molecule]['xlogP']
		if xlogP == None or xlogP == '':
			return None

		scaled_xlogP = (xlogP-self.mapDict['xlogP']['min'])/(self.mapDict['xlogP']['max']-self.mapDict['xlogP']['min'])
		return scaled_xlogP

	def getPKA(self):
		pass

	def getSolubility(self):
		pass

	def getSize(self):
		pass

	def getStability(self):
		pass

	def getWeight(self):
		'''molecular_weight = self.molecule.molecular_weight
		if molecular_weight == None:
			self.cannotUseInPredictiveModel = True
			return None
		molecular_weight = float(molecular_weight)'''

		if self.molecule not in self.drugsDict:
			return None

		molecular_weight = self.drugsDict[self.molecule]['molecular_weight']

		if molecular_weight ==  None or molecular_weight == '':
			return None

		molecular_weight_scaled = (molecular_weight-self.mapDict['molecular_weight']['min'])/(self.mapDict['molecular_weight']['max'] - self.mapDict['molecular_weight']['min'])
		return molecular_weight_scaled

	def getMass(self):
		'''mass = self.molecule.exact_mass
		if mass == None:
			self.cannotUseInPredictiveModel = True
			return None
		mass = float(mass)'''
		if self.molecule not in self.drugsDict:
			return None

		mass = self.drugsDict[self.molecule]['mass']

		if mass == None or mass == '':
			return None
		mass_scaled = (mass-self.mapDict['mass']['min'])/(self.mapDict['mass']['max']-self.mapDict['mass']['min'])
		return mass_scaled

	#Convert compound name to smiles
	#Note that this needs some attention. This is such clean code.
	def getSmiles(self):
		'''if self.molecule == None:
			self.cannotUseInPredictiveModel = True
			return None
		self.smiles = self.molecule.canonical_smiles'''

		if self.molecule not in self.drugsDict:
			return None

		self.smiles = self.drugsDict[self.molecule]['SMILES']

		if self.smiles == None or self.smiles == '':
			return None

		return self.smiles


	def molToVec(self):
		'''if self.molecule == None:
			return None
		if self.smiles == None:
			self.getSmiles()'''

		if self.smiles == None:
			self.getSmiles()
			if self.smiles == None:
				return None

		molMol = Chem.MolFromSmiles(self.smiles)
		sentence = MolSentence(mol2alt_sentence(molMol, radius=1))
		self.vector = sentences2vec([sentence], self.model, unseen='UNK')
		self.vector = DfVec(self.vector[0]).vec
		return self.vector






