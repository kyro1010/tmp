+import os
import boto3
import botocore
import pandas as pd
import time
import mysql.connector

import numpy as np
from compoundPrep import MoleculePrepper
import random
from pharmaDb import engine, maxVecNormTbl
import math

conn = engine.connect()


def scalarNoise(scalarEpsilon):
	return np.random.normal(0,scalarEpsilon)



def smilesToVec(smiles):
	pass


def vecToStr(vec):
	vecStr = ""
	for elmt in vec:
		vecStr += str(elmt)
		vecStr += ","

	return vecStr


def strToVec(vecStr):
	vec = []
	vecStrList = vecStr.split(',')
	for elmt in vecStrList:
		vec.append(float(elmt))

	return vec


def diffVec(v1,v2):
	diff = []
	for i in range(len(v1)):
		diff.append(0)
	for i in range(len(v1)):
		diff[i] = v1[i] - v2[i]

	return diff


def vecNorm(v):
	norm = 0
	for e in v:
		norm += float((e*e))
		#print(e*e)

	return float(norm**0.5)

def printVec(v):
	vecStr = ""
	for e in v:
		vecStr += (str(e) + ", ")
	#print(vecStr)

def addNoiseToVec(seedVector, radius, stepSize):
	resultantVector = [x for x in seedVector]	#We will modify resultantVector.

	numDimensions = len(seedVector)

	indexInfo = dict()
	#vecN = vecNorm(diffVec(resultantVector,seedVector))
	vecN = 0
	while vecN < radius:
		#Select index at random.
		randIndex = random.randint(0,numDimensions-1)
		#print(randIndex)
		randSign = None

		if randIndex not in indexInfo:
			randSign = random.randint(0,1)
			indexInfo[randIndex] = randSign
		else:
			randSign = indexInfo[randIndex]

		if randSign == 0:
			resultantVector[randIndex] = resultantVector[randIndex] - stepSize
		else:
			resultantVector[randIndex] = resultantVector[randIndex] + stepSize


		vecN = vecNorm(diffVec(resultantVector,seedVector))
		#printVec(diffVec(resultantVector, seedVector))
	return resultantVector


def vecMagnitude(vec):
	if vec == None or len(vec) == 0:
		return None

	ss = 0

	for elmt in vec:
		ss += (elmt**2)

	return ss**0.5






while True:
        #Set up session.
        session = boto3.session.Session()
        client = session.client('s3',
                        endpoint_url='https://nyc3.digitaloceanspaces.com', # Find your endpoint in the control p>
                        config=botocore.config.Config(s3={'addressing_style': 'virtual'}), # Configures to use su>
                        region_name='nyc3', # Use the region in your endpoint.
                        aws_access_key_id='DO00V9QUUNT6EYJFW8LP', # Access key pair. You can create access key pa>
                        aws_secret_access_key='8xM5IsfhK8GJY47hWdFB0/7Q7UDXksTM1K/v3kX3px4') # Secret access key >


        #Query the database. In particular the 'drugsOriginal' table.
        mydb = mysql.connector.connect(
        host ="dbaas-db-6613674-do-user-13732574-0.b.db.ondigitalocean.com",
        user = "doadmin",
        password ="AVNS_HYWrZBGRK742SWOrp31",
        port = 25060,
        database = "test3"
        )

        mycursor = mydb.cursor()

	'''We need to augment the drug data. We need to get the parameters from the dataAugmentation table. There are parameters which specify the number of synthetic datum points for each drug,
	as well as the epsilon window for the molecular vectors as well as the feature scalars.'''

	mycursor.execute("SELECT * FROM augmentedDataParameters")
	rows = mycursor.fetchall()
	if len(rows) == 0:
		mycursor.execute("INSERT INTO augmentedDataParameters (numAugmentedDatumPoints, vecEpsilon, scalarEpsilon) VALUES ('%i','%f','%f')"%(DEFAULT_NUM_AUGMENTED_DATUM_POINTS, DEFAULT_VEC_EPSILON, DEFAULT_SCALAR_EPSILON))
		mydb.commit()

	mycursor.execute("SELECT * FROM augmentedDataParameters")
	row = mycursor.fetchall()[0]
	numAugmentedDatumPoints = row[1]
	vecEpsilon = row[2]
	scalarEpsilon = row[3]

	#Select all drugs from the synthetic drugs table.
	mycursor.execute("SELECT * FROM drugsSynthetic")
	rows = mycursor.fetchall()

	augmentedParentDrugsSet = set()
	for row in rows:
		parentDrugId = row[1]
		augmentedParentDrugsSet.add(parentDrugId)

	#Now iterate over all parent drugs and for any that has not yet had an augmented set generated for it, produce the augmented set.
	mycursor.execute("SELECT * FROM drugsOriginal")
	rows = mycursor.fetchall()

	for row in rows:
		parentDrugId = row[0]
		molVector = row[4]
		absorption = row[5]
		distribution = row[6]
		metabolism = row[7]
		excretion = row[8]
		toxicity = row[9]
		xlogp = row[10]
		pka = row[11]
		solubility = row[12]
		size = row[13]
		stability = row[14]
		weight = row[15]
		mass = row[16]

		molVector = strToVec(molVector)

		parentScalarFeaturesList = [absorption, distribution, metabolism, excretion, toxicity, xlogp, pka, solubility, size, stability, weight, mass]
		scalarFeatureNames = ['absorption', 'distribution', 'metabolism', 'excretion', 'toxicity', 'xlogp', 'pka', 'solubility', 'size', 'stability', 'weight', 'mass']

		if parentDrugId not in augmentedParentDrugsSet:
			#Must perform data augmentation.
			for i in range(numAugmentedDatumPoints):
				#Produce augmented vector.
				scalarFeaturesDict = dict()
				scalarFeaturesDict['absorption'] = None
				scalarFeaturesDict['distribution'] = None
				scalarFeaturesDict['metabolism'] = None
				scalarFeaturesDict['excretion'] = None
				scalarFeaturesDict['toxicity'] = None
				scalarFeaturesDict['xlogp'] = None
				scalarFeaturesDict['pka'] = None
				scalarFeaturesDict['solubility'] = None
				scalarFeaturesDict['size'] = None
				scalarFeaturesDict['stability'] = None
				scalarFeaturesDict['weight'] = None
				scalarFeaturesDict['mass'] = None

				#Produce augmented scalars.
				for feature, featureName in zip(parentScalarFeaturesList, scalarFeatureNames):
					#TODO Produce augmentations here.
					if feature == None:
						continue

					scalarFeaturesDict[featureName] = feature + scalarNoise(scalarEpsilon)

				#Produce augmented vectors.
				molVectorPrime = addNoiseToVec(molVec, vecEpsilon, 0.2)
				molVectorPrimeStr = vecToStr(molVectorPrime)

				#TODO Store this data in the database.
				mycursor.execute("INSERT INTO drugsSynthetic (parentId, vector, absorption, distribution, metabolism, excretion, toxicity, xlogP, pKa, solubility, size, stability, molecularWeight, mass) VALUES ('%i', '%s', '%f','%f', '%f','%f','%f','%f', '%f','%f','%f','%f','%f','%f')"%(parentDrugId, molVectorPrimeStr, scalarFeaturesDict['absorption'],scalarFeaturesDict['distribution'],scalarFeaturesDict['metabolism'],scalarFeaturesDict['excretion'], scalarFeaturesDict['toxicity'],scalarFeaturesDict['xlogp'], scalarFeaturesDict['pka'], scalarFeaturesDict['solubility'], scalarFeaturesDict['size'], scalarFeaturesDict['stability'], scalarFeaturesDict['weight'], scalarFeaturesDict['mass']))
				mydb.commit()







