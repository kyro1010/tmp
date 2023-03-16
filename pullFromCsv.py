import os
import boto3
import botocore
import pandas as pd
import time
import mysql.connector
import numpy as np
import random
import math

import pubchempy as pc
from rdkit import Chem
from mol2vec.features import mol2alt_sentence
import pybel
from gensim.models import word2vec
from mol2vec.features import mol2alt_sentence, mol2sentence, MolSentence, DfVec, sentences2vec

#from gensim.models import Word2Vec

#model = Word2Vec.load('model_300dim.pkl')
#model = word2vec.Word2Vec.load('mol2vec/examples/models/model_300dim.pkl')
model = word2vec.Word2Vec.load('model_300dim.pkl')

def molToVec(smiles):
	if smiles == None:
		return None

	molMol = Chem.MolFromSmiles(smiles)
	sentence = MolSentence(mol2alt_sentence(molMol, radius=1))
	vector = sentences2vec([sentence], model, unseen='UNK')
	vector = DfVec(vector[0]).vec
	return vector


def vecToStr(vec):
	vecStr = ""
	for elmt in vec:
		vecStr += str(elmt)
		vecStr += ","

	return vecStr

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

	#mycursor.execute("SELECT * FROM drugsOriginal")
	#rows = mycursor.fetchall()
	#drugDict = dict()



	#Set the name of your Space and the file you want to download
	space_name = 'automatedchemist'
	file_name = 'drugsData.csv'

	# Download the file from your Space
	client.download_file(space_name, file_name, file_name)

	#TODO Add in some error handling here.
	#Pull csv file from the bucket and store on the disk.
	df = pd.read_csv(file_name)
	for index, row in df.iterrows():
		disease = row['Disease'].lower().strip()
		target = row['Target'].lower().strip()
		drug = row['Drug'].lower().strip()
		efficacy = float(row['Efficacy'])
		print(disease)

		#Check to see if the drug is in the database.
		mycursor.execute("SELECT * FROM drugsOriginal WHERE drugsOriginal.drugName = '%s'"%(drug))
		rows = mycursor.fetchall()
		if len(rows) == 0:
			#Insert this data in the database.
			#TODO: Pull features data from pubchem.
			smiles = None
			xlogP = None
			molecularWeight = None
			absorption = None
			distribution = None
			metabolism = None
			excretion = None
			toxicity = None
			xlogP = None
			pKa = None
			solubility = None
			size = None
			stability = None
			molecularWeight = None
			mass = None

			molVector = ""

			pubChemDrugIds = pc.get_cids(drug)

			mySqlFeatures = ['drugName'] #To be built upon.
			mySqlFeaturesActual = [drug] #To be built upon.
			stringFormattingList = ['%s'] #To be built upon.


			if len(pubChemDrugIds) > 0:
				pubChemDrugId = pubChemDrugIds[0]
				d = pc.get_properties(['canonical_smiles', 'xlogp','exact_mass','molecular_weight'], pubChemDrugId)
				if len(d) > 0:
					d = d[0]
					if 'CanonicalSMILES' in d:
						smiles = d['CanonicalSMILES']
						mySqlFeatures.append('SMILES')
						mySqlFeaturesActual.append(smiles)
						stringFormattingList.append('%s')
						#drugInfo[id_]['SMILES'] = None

					if 'XLogP' in d:
						xlogP = d['XLogP']
						mySqlFeatures.append('base_xlogP')
						mySqlFeaturesActual.append(xlogP)
						stringFormattingList.append('%f')
						#drugInfo[id_]['baseXLogP'] = None

					if 'ExactMass' in d:
						mass = float(d['ExactMass'])
						mySqlFeatures.append('base_mass')
						mySqlFeaturesActual.append(mass)
						stringFormattingList.append('%f')
						#drugInfo[id_]['baseMass'] = None

					if 'MolecularWeight' in d:
						molecularWeight = float(d['MolecularWeight'])
						mySqlFeatures.append('base_molecularWeight')
						mySqlFeaturesActual.append(molecularWeight)
						stringFormattingList.append('%f')
						#drugInfo[id_]['baseMolecularWeight'] = None


					#TODO: Then vectorize the molecule from the SMILES representation.
					drugVector = molToVec(smiles)
					#TODO: Now  convert that vector to a string.
					molVector = vecToStr(molVector)

					mySqlFeatures.append('vector')
					mySqlFeaturesActual.append(molVector)
					stringFormattingList.append('%s')

			mySqlFeaturesString = "("
			for i in range(len(mySqlFeatures)-1):
				mySqlFeaturesString += mySqlFeatures[i]
				mySqlFeaturesString += ","
			mySqlFeaturesString += mySqlFeatures[-1]
			mySqlFeaturesString += ")"

			stringFormattingString = "("
			for i in range(len(stringFormattingList)-1):
				stringFormattingString += "\'"+stringFormattingList[i]+"\'"
				stringFormattingString += ","
			stringFormattingString += "\'"+stringFormattingList[-1]+"\'"
			stringFormattingString += ")"

			mySqlFeaturesActual = tuple(mySqlFeaturesActual)
			print(stringFormattingString)
			print(("INSERT INTO drugsOriginal "+ mySqlFeaturesString+" VALUES "+ str(stringFormattingString)) % mySqlFeaturesActual)
			mycursor.execute(("INSERT INTO drugsOriginal "+ mySqlFeaturesString+" VALUES "+ str(stringFormattingString)) % mySqlFeaturesActual)
			mydb.commit()

		#Check to see if the disease is in the database.
		mycursor.execute("SELECT * FROM diseases WHERE diseases.diseaseName = '%s'"%(disease))
		rows = mycursor.fetchall()
		if len(rows) == 0:
			#Insert this data into the database.
			mycursor.execute("INSERT INTO diseases (diseaseName) VALUES ('%s')"%(disease))
			mydb.commit()

		#Check to see if the target is in the database.
		mycursor.execute("SELECT * from targets WHERE targets.targetName = ('%s')"%(target))
		rows = mycursor.fetchall()
		if len(rows) == 0:
			#Insert this data into the database.
			mycursor.execute("INSERT INTO targets (targetName) VALUES ('%s')"%(target))
			mydb.commit()

		#Get ids.
		mycursor.execute("SELECT * FROM drugsOriginal WHERE drugsOriginal.drugName = '%s'"%(drug))
		row = mycursor.fetchall()[0]
		drugId = row[0]

		mycursor.execute("SELECT * FROM diseases WHERE diseases.diseaseName = '%s'"%(disease))
		row = mycursor.fetchall()[0]
		diseaseId = row[0]

		mycursor.execute("SELECT * from targets WHERE targets.targetName = '%s'"%(target))
		row = mycursor.fetchall()[0]
		targetId = row[0]

		mycursor.execute("SELECT * from drugDiseaseTargetTuplesOriginal WHERE drugDiseaseTargetTuplesOriginal.parentDrugId = ('%i') and drugDiseaseTargetTuplesOriginal.targetId = ('%i') and drugDiseaseTargetTuplesOriginal.diseaseId = ('%i')"%(drugId,targetId,diseaseId))
		rows = mycursor.fetchall()
		#First check to see if the drug is in the database. Then check to see if the target is, then check to see if the disease is.
		if len(rows) == 0:
			#Insert into the database.
			mycursor.execute("INSERT INTO drugDiseaseTargetTuplesOriginal (diseaseId, targetId, parentDrugId, efficacy) VALUES ('%i','%i','%i','%f')"%(diseaseId,targetId,drugId,efficacy))
			mydb.commit()



		#Commit these changes.
		mydb.commit()


	#Extract data from the file. Determine if any new data has been added. If so, store the new data in the database.

	time.sleep(60)
