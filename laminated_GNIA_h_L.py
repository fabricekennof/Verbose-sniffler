from abaqus import *
from abaqusConstants import *
from caeModules import *
from visualization import XYData,USER_DEFINED
from odbAccess import *
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import csv
import time
import math
import subprocess
from collections import defaultdict
import os
import methodCallback
from job import JobType
import re
import fnmatch

#transferrable units in MM-MPA-etc..


Analysis_type='GNA_h_L' #Define general name for analysis type
dir="C:\\Temp_fakennof\%s"%Analysis_type #maak een map aan

n=1
while (os.path.isdir(dir)):
	Analysis_type='GNA_h_L_%s' %str(n)
	dir="C:\\Temp_fakennof\%s"%Analysis_type
	n+=1

os.mkdir(dir)
os.chdir(r"%s"%dir)

Beam=dict()
Laminate=dict()

###########################################################################
#_____Analysis Choices____################################################################
###########################################################################

# 0 is present, 1  is not.
N_beams_list=[2,3,4] #N_beams is defined as a beam with N equal to the amount of glass beams 
lateral_restraint=0 
Viscoelastic=0
Load_temperature=20
interlayer_material='SGP'#PVB or SGP

#Type glass
Glass_type='ANG' #ANG or HSG of FTG

#General remark: in XZ plane, all the elements should have equal nodes
Mesh_size=20 #Define the necessary sizes of the mesh in mm
Beam_element=C3D20R #Define the elements for the glass panes, note that for C3D8R hourglassControl=ENHANCED
Interlayer_element=C3D20 #Up to present interlayer had to be defined, even though it is not used for monolithic beams

#Dimensions
Interlayer_thickness=1.52
Beam_thickness=15
L_hs=[8,10,12,15,20,25] #weggelaten voorlopig
Beam_length_total=[5000,6000]#1000,2000 tijdelijk verwijderd
Beam_lengths=map(lambda x:x*0.5,Beam_length_total)#due to symmetry



#Lateral restraints
##k_sil=1 #define spring contant for lateral restraint N/mm/mm
##K_sil=k_sil*Beam_length

#Imperfections
impSizes=[200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]     


###########################################################################
#_____partial scripts____################################################################
###########################################################################

def LBA(name):
    nameStep='LTB'
    nameOdb=name+'.odb'
    odb=session.openOdb(name=nameOdb)
    step=odb.steps[nameStep]
    descr=step.frames[1].description
    eigenval=float(descr.split('=')[-1])
    return eigenval
    
def Calculation_time(name):
    m,s =divmod(name,60)
    h,m=divmod(m,60)
    return "%02d:%02d:%02d" % (h, m, s)
    
def amount_elements(parameter):
    parameter=int(parameter)
    amount_beam=0
    amount_interlayer=0
    total_amount=0

    amount_beam=len((mdb.models[file_name].parts['Beam'].elements))*(parameter)
    if(parameter>1):
        amount_interlayer=len(mdb.models[file_name].parts['Interlayer'].elements)*(parameter)-len(mdb.models[file_name].parts['Interlayer'].elements)
   
    
    total_amount=amount_beam+amount_interlayer
    return total_amount
    
def addLinesAfter(fileNameBron,fileNameRes,lines,afterLine):
	fileBron=open(fileNameBron,'r')
	fileRes=open(fileNameRes,'w+')
	for line in fileBron:
		if line.find(afterLine)>=0:
			fileRes.write(line)
			fileRes.write(lines)
		else:
			fileRes.write(line)	
	fileBron.close()
	fileRes.close()


def copyContent(sourceFileName,destinationFile):
	# kopieer de inhoud van een bestand sourcefileName naar een reeds geopende destinationFile
	sourceFile=open(sourceFileName,'r')
	lines = sourceFile.readlines()
	for line in lines:
		destinationFile.write(line)
	sourceFile.close()
	
def Glass_stress(Glass_type):
	if (Glass_type=='ANG'):
		return 45
	elif (Glass_type=='HSG'):
		return 70
	elif (Glass_type=='FTG'):
		return 120

def F_cr(N_beams,Beam_thickness,Beam_height,Beam_length):
    Beam_length=Beam_length*2
    Beam_thickness=N_beams*Beam_thickness
    I=Beam_height*math.pow(Beam_thickness,3)/12.
    I_t=I*4
    E=70000
    G=28456
    c1=1.365
    c2=0.553
    z_a=Beam_height/2.
    
    F_cr=c1*4*math.pow(math.pi,2)*E*I/math.pow(Beam_length,3)*(math.sqrt(math.pow(Beam_length,2)*G*I_t/(math.pow(math.pi,2)*E*I)+pow(c2*z_a,2))-c2*z_a)
    return F_cr
    
def Estimate_location(node_number,instance_name,Beam_length,Beam_height):
	x=odb.rootAssembly.instances[instance_name].nodes[node_number].coordinates[0]
	y=odb.rootAssembly.instances[instance_name].nodes[node_number].coordinates[1]
	
	if(0<=x and x<Beam_length/4.):
		locatie_x='Centre'
	elif(x==Beam_length/4.):
		locatie='Warning, edge of Centrezone'
	elif(Beam_length/4.<x and x<3*Beam_length/4.):
		locatie_x='Zone between support and centre'
	elif(x==3*Beam_length/4.):
		locatie_x='Warning, edge of Support zone'
	elif(3*Beam_length/4.<x and x<=Beam_length):
		locatie_x='Support zone'
	
	if(0<=y and y<Beam_height/3.):
		locatie_y='Bottom'
	elif(y==Beam_height/3.):
		locatie_y='Warning: border between bottom and middle'
	elif(Beam_height/3.<y and y<2*Beam_height/3.):
		locatie_y='Middle'
	elif(y==2*Beam_height/3.):
		locatie_y='Warning: border between top and middle'
	elif(2*Beam_height/3.<y and y<=Beam_height):
		locatie_y='Top'
	locatie_instance='%s'%instance_name
	
	return ('x_%s,y_%s,inst_%s')%(locatie_x,locatie_y,locatie_instance)

def newestStat():
	lijst_stats=[]
	for file in os.listdir('.'):
		if fnmatch.fnmatch(file, '*.sta'):
			modtime=os.stat(file).st_mtime
			lijst_stats.append([file,modtime])
	if len(lijst_stats)>0:
		lijst_stats2=sorted(lijst_stats, key=lambda x:x[1])
		sta=lijst_stats2[-1][0]
	else:
		sta='none.sta'
	return sta

def ODB_open_safe(odb):
	for x in range(0, 10):  # try 10 times
		try:
			odb=session.openOdb(name='%s.odb'%odb)
			str_error = None
		except OdbError,e:
			str_error=str(e)
			pass
		except IOError,e:
			str_error=str(e)
			pass
		if str_error:
			print str_error
			time.sleep(2)  # wait for 2 seconds before trying to fetch the data again
		else:
			return odb
			break
			

def send_email(user, pwd, recipient, subject, body):
    import smtplib

    gmail_user = user
    gmail_pwd = pwd
    FROM = user
    TO = recipient if type(recipient) is list else [recipient]
    SUBJECT = subject
    TEXT = body

    # Prepare actual message
    message = """From: %s\nTo: %s\nSubject: %s\n\n%s
    """ % (FROM, ", ".join(TO), SUBJECT, TEXT)
    try:
        server = smtplib.SMTP("smtp.gmail.com", 587)
        server.ehlo()
        server.starttls()
        server.login(gmail_user, gmail_pwd)
        server.sendmail(FROM, TO, message)
        server.close()
        print 'successfully sent the mail'
    except:
        print "failed to send mail"
    
###########################################################################
#_____model____#####################################################################
###########################################################################

Plotnr=1 #for plotting XY-data properly
Chartnr=1

file=open(Analysis_type+'.csv','ab')
fileResults=csv.writer(file)
fileResults.writerow(['file_name,Calctime,amount_of_elements, Buckling_resistance_ANG,Endframe_ANG,Breakstress_ANG,Fail_instance_ANG,Fail_element_ANG, Estimate_location_ANG,Buckling_resistance_HSG,Endframe_HSG,Breakstress_HSG,Fail_instance_HSG,Fail_element_HSG,Estimate_location_HSG,Buckling_resistance_FTG,Endframe_FTG,Breakstress_FTG,Fail_instance_FTG,Fail_element_FTG,Estimate_location_FTG'])
file.close()

for N_beams in N_beams_list:
	for Beam_length in Beam_lengths:
		for L_h in L_hs:
			for impSize in impSizes:
				Total_thickness=N_beams*Beam_thickness+(N_beams-1)*Interlayer_thickness
				file_name=Analysis_type+"_Beams_%s_Beam_length_%s_Lh_%s_imp_%s" %(int(N_beams),int(Beam_length), int(L_h),int(impSize))
				print file_name
				Beam_height=Beam_length*2.0/L_h
				imp=Beam_length*2.0/impSize
				
				instances=[] #een lijst aanmaken met de instances voor de imperfecties
				instances_glass=[]#Voor spanningscontrole
				instances_interlayer=[]
			
				dir_temp="%s\%s"%(dir,file_name)
				os.mkdir(dir_temp)
				os.chdir(r"%s"%dir_temp)
				
				##################
				#Create MODEL DATA BASE#
				##################
				newmodel=Mdb()
				mdb.Model(name=file_name)
				mdb.models[file_name].setValues(absoluteZero=0)
				del mdb.models['Model-1']
				
				#_______________Create PART______________#
				
				#Sketch
				s1 = mdb.models[file_name].ConstrainedSketch(name='__profile__', sheetSize=200.0)
				g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
				s1.setPrimaryObject(option=STANDALONE)
				s1.rectangle(point1=(0, 0), point2=(Beam_length, Beam_height))
				
				Beam= mdb.models[file_name].Part(name='Beam', dimensionality=THREE_D, type=DEFORMABLE_BODY)
				Beam.BaseSolidExtrude(sketch=s1, depth=Beam_thickness)

				Interlayer = mdb.models[file_name].Part(name='Interlayer', dimensionality=THREE_D, type=DEFORMABLE_BODY)
				Interlayer.BaseSolidExtrude(sketch=s1, depth=Interlayer_thickness)
				
				#Material assignment
				mdb.models[file_name].Material(name='Glass')
				mdb.models[file_name].materials['Glass'].Elastic(table=((70000.0, 0.23), ))
				mdb.models[file_name].HomogeneousSolidSection(name='Beamsection', material='Glass',)
			    
				if(Viscoelastic==0):
					if(interlayer_material=='PVB'):
						mdb.models[file_name].Material(name='PVB')
						mdb.models[file_name].materials['PVB'].Elastic(table=((18, 0.45), )) #so high because buckling loads are usually quite short in time, temperature not specifically taken into account
						mdb.models[file_name].HomogeneousSolidSection(name='Interlayersection', material='PVB',)
					elif(interlayer_material=='SGP'):
						mdb.models[file_name].Material(name='SGP')
						mdb.models[file_name].materials['SGP'].Elastic(table=((149, 0.48), ))
						mdb.models[file_name].HomogeneousSolidSection(name='Interlayersection', material='SGP',)
					
				elif(Viscoelastic==1):
					if(interlayer_material=='PVB'):
						mdb.models[file_name].Material(name='PVB')
						mdb.models[file_name].materials['PVB'].Elastic(moduli=INSTANTANEOUS,table=((471, 0.45), ))
						mdb.models[file_name].materials['PVB'].Viscoelastic(domain=TIME,frequency=PRONY, table=((0.1606, 0.0, 3.2557e-11), (0.07877, 0.0, 4.9491e-09), (0.2912, 0.0, 7.2427e-08), (0.071155, 0.0, 9.8635e-06), (0.2688, 0.0, 0.0028059), (0.089568, 0.0, 0.16441), (0.030183, 0.0, 2.2648), (0.0076056, 0.0, 35.364), (0.0009634, 0.0, 9367.5), (0.0004059, 0.0, 641410.0), (0.0006143, 0.0, 41347000.0))) #g_i,k_i_tau_i
						mdb.models[file_name].materials['PVB'].viscoelastic.Trs(table=((20.0, 20.7, 91.1), )) #Theta_0,C1,C2
						mdb.models[file_name].HomogeneousSolidSection(name='Interlayersection', material='PVB',)
					elif(interlayer_material=='SGP'):
						mdb.models[file_name].Material(name='SGP')
						mdb.models[file_name].materials['SGP'].Elastic(moduli=INSTANTANEOUS,table=((447, 0.49), ))
						mdb.models[file_name].materials['SGP'].Viscoelastic(domain=TIME,frequency=PRONY, table=((0.5932, 0.0, 0.065173), (0.1122, 0.0, 0.9669), (-0.0051988, 0.0, 82.31), (0.049333, 0.0, 446.3), (0.020831, 0.0, 5648.0), (0.061392, 0.0, 65132.0), (0.043697, 0.0, 504060.0), (0.050251, 0.0, 4908400.0), (0.029005, 0.0, 3.3452e-07), (0.019283, 0.0, 523630000.0), (0.007369, 0.0, 7739600000.0), (0.0054495, 0.0, 126130000000.0), (0.013155, 0.0, 8331600000000.0))) #g_i,k_i_tau_i
						mdb.models[file_name].materials['SGP'].viscoelastic.Trs(table=((20.,135., 760.), )) #Theta_0,C1,C2
						mdb.models[file_name].HomogeneousSolidSection(name='Interlayersection', material='SGP',)
		    
				#Assign section
				Cell1=Interlayer.cells.getByBoundingBox(xMin=0,xMax=Beam_length,yMin=0,yMax=Beam_height)
				Interlayer.Set(name='Interlayerset',cells=Cell1)
				Interlayer.SectionAssignment(region=Interlayer.sets['Interlayerset'],sectionName='Interlayersection', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
				
				Cell2=Beam.cells.getByBoundingBox(xMin=0,xMax=Beam_length,yMin=0,yMax=Beam_height)
				Beam.Set(name='Beamset',cells=Cell2)
				Beam.SectionAssignment(region=Beam.sets['Beamset'],sectionName='Beamsection', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
				
				#Mesh
				if (Beam_element==C3D8R):
					elemBeam = mesh.ElemType(elemCode=Beam_element, elemLibrary=STANDARD, hourglassControl=ENHANCED)
				else:
					elemBeam = mesh.ElemType(elemCode=Beam_element, elemLibrary=STANDARD)
				elemInterlayer=mesh.ElemType(elemCode=Interlayer_element, elemLibrary=STANDARD)
				    
				Mesh_length_division=math.ceil(Beam_length/(Mesh_size*1.0))#gives more control over mesh generation
				Mesh_height_division=math.ceil(Beam_height/(Mesh_size*1.0))
				    
				Seed_length_beam=mdb.models[file_name].parts['Beam'].edges.findAt(((Beam_length/2.,0,0),))
				Seed_height_beam=mdb.models[file_name].parts['Beam'].edges.findAt(((0,Beam_height/2.,0),))
				mdb.models[file_name].parts['Beam'].setElementType(regions=mdb.models[file_name].parts['Beam'].sets['Beamset'], elemTypes=(elemBeam,))
				mdb.models[file_name].parts['Beam'].seedEdgeByNumber(edges=Seed_length_beam, number=int(Mesh_length_division), constraint=FINER)
				mdb.models[file_name].parts['Beam'].seedEdgeByNumber(edges=Seed_height_beam, number=int(Mesh_height_division), constraint=FINER)
				mdb.models[file_name].parts['Beam'].seedPart(deviationFactor=0.1, size=Mesh_size)
				mdb.models[file_name].parts['Beam'].generateMesh()
				    
				Seed_thickness_interlayer=mdb.models[file_name].parts['Interlayer'].edges.findAt(((0,0,Interlayer_thickness/2.),))
				Seed_length_interlayer=mdb.models[file_name].parts['Interlayer'].edges.findAt(((Beam_length/2.,0,0),))
				Seed_height_interlayer=mdb.models[file_name].parts['Interlayer'].edges.findAt(((0,Beam_height/2.,0),))
				mdb.models[file_name].parts['Interlayer'].setElementType(regions=mdb.models[file_name].parts['Interlayer'].sets['Interlayerset'], elemTypes=(elemInterlayer,))
				mdb.models[file_name].parts['Interlayer'].seedEdgeByNumber(edges=Seed_thickness_interlayer, number=2, constraint=FINER)
				mdb.models[file_name].parts['Interlayer'].seedEdgeByNumber(edges=Seed_length_interlayer, number=int(Mesh_length_division), constraint=FINER)
				mdb.models[file_name].parts['Interlayer'].seedEdgeByNumber(edges=Seed_height_interlayer, number=int(Mesh_height_division), constraint=FINER)
				mdb.models[file_name].parts['Interlayer'].seedPart(deviationFactor=0.1, size=Mesh_size, minSizeFactor=0.1)
				mdb.models[file_name].parts['Interlayer'].generateMesh()
				    
				exclusion_edge=mdb.models[file_name].parts['Beam'].edges.findAt(((0,Beam_height,Beam_thickness/2.),),((Beam_length,0,Beam_thickness/2.),))
					    
				exclusion_set=mdb.models[file_name].parts['Beam'].Set(edges=exclusion_edge, name='exclusion_set')
				Stress_exclusion_elements=exclusion_set.elements
				Stress_exclusion=mdb.models[file_name].parts['Beam'].Set(elements=Stress_exclusion_elements,name='Stress_exclusion')
				All_elements=elements=mdb.models[file_name].parts['Beam'].elements.getByBoundingBox(xMin=0,xMax=Beam_length,yMin=0,yMax=Beam_height,zMin=0,zMax=Beam_thickness)    
				All_elements=mdb.models[file_name].parts['Beam'].Set(elements=All_elements,name='All_elements')
				Stress_inclusion=mdb.models[file_name].parts['Beam'].SetByBoolean(name='Stress_inclusion', operation=DIFFERENCE, sets=(All_elements,Stress_exclusion,))
				
				
				#_______________Create Instance and Necessary interactions______________#
				mdb.models[file_name].rootAssembly.DatumCsysByDefault(CARTESIAN)
				
				Face_links_interlayer=dict()
				Face_links_beam=dict()
				
				Face_rechts_interlayer=dict()
				Face_rechts_beam=dict()
				Translation_beam=0
				Translation_interlayer=0
				
				for x in range (1,N_beams+1):
					Translation_beam=(x-1)*(Interlayer_thickness+Beam_thickness)
					mdb.models[file_name].rootAssembly.Instance(name='Plate_%s'%x, part=Beam, dependent=ON)
					instances.append(mdb.models[file_name].rootAssembly.instances['Plate_%s'%x])
					instances_glass.append(mdb.models[file_name].rootAssembly.instances['Plate_%s'%x])
					mdb.models[file_name].rootAssembly.translate(instanceList=('Plate_%s'%x, ), vector=(0.0, 0.0, Translation_beam))
				    
					Face_links_beam[x]=mdb.models[file_name].rootAssembly.instances['Plate_%s'%x].faces.getByBoundingBox(xMin=0,xMax=Beam_length,yMin=0,yMax=Beam_height,zMin=Translation_beam,zMax=Translation_beam)
					Face_rechts_beam[x]=mdb.models[file_name].rootAssembly.instances['Plate_%s'%x].faces.getByBoundingBox(xMin=0,xMax=Beam_length,yMin=0,yMax=Beam_height,zMin=Translation_beam+Beam_thickness,zMax=Translation_beam+Beam_thickness)
				
					if(x>=2):
						Translation_interlayer=(x-2)*(Interlayer_thickness+Beam_thickness)+Beam_thickness
						mdb.models[file_name].rootAssembly.Instance(name='Interlayer_%s'%(x-1), part=Interlayer, dependent=ON)
						mdb.models[file_name].rootAssembly.translate(instanceList=('Interlayer_%s'%(x-1), ), vector=(0., 0., Translation_interlayer))
						instances.append(mdb.models[file_name].rootAssembly.instances['Interlayer_%s'%(x-1)])                  
						Face_links_interlayer[x-1]=mdb.models[file_name].rootAssembly.instances['Interlayer_%s'%(x-1)].faces.getByBoundingBox(xMin=0,xMax=Beam_length,yMin=0,yMax=Beam_height,zMin=Translation_interlayer,zMax=Translation_interlayer)
						Face_rechts_interlayer[x-1]=mdb.models[file_name].rootAssembly.instances['Interlayer_%s'%(x-1)].faces.getByBoundingBox(xMin=0,xMax=Beam_length,yMin=0,yMax=Beam_height,zMin=Translation_interlayer+Interlayer_thickness,zMax=Translation_interlayer+Interlayer_thickness)
					
						mdb.models[file_name].Tie(name='Constraint_%s'%(2*x-3), master=mdb.models[file_name].rootAssembly.Surface(side1Faces=Face_rechts_beam[x-1], name='Plate_%s_region_right'%(x-1)), slave=mdb.models[file_name].rootAssembly.Surface(side1Faces=Face_links_interlayer[x-1], name='Interlayer_%s_links_region'%(x-1)), positionToleranceMethod=COMPUTED, adjust=OFF, tieRotations=ON, thickness=ON)
						mdb.models[file_name].Tie(name='Constraint_%s'%(2*x-2), master=mdb.models[file_name].rootAssembly.Surface(side1Faces=Face_links_beam[x], name='Plate_%s_region_left'%(x)), slave=mdb.models[file_name].rootAssembly.Surface(side1Faces=Face_rechts_interlayer[x-1], name='Interlayer_%s_rechts_region'%(x-1)), positionToleranceMethod=COMPUTED, adjust=OFF, tieRotations=ON, thickness=ON)
				    
				#____________Imperfection introduction_______________#   
				
				
				impNaam=file_name.split('_')[0]+'_'+file_name.split('_')[1]+'_'+file_name.split('_')[2]+'_'+file_name.split('_')[3]+'_'+file_name.split('_')[4]+'_imp'+str(int(impSize))
			    
				inpImp=open(impNaam+'.inp','w')
				for i in instances:
					for node in i.nodes:
						x=node.coordinates[0]
						y=node.coordinates[1]
						vz=imp*(math.cos(math.pi*x/(Beam_length*2.)))
						label=node.label
						inpImp.write('%s.%s,0., 0., %s\n' %(i.name,label,vz))
				inpImp.close()
		    
				#____________Boundary conditions_______________# 
				Beam_ur1=mdb.models[file_name].rootAssembly.instances['Plate_1'].edges.findAt(((Beam_length,Beam_height/float(2),0.),))
				edges_RBC_ur1=Beam_ur1+mdb.models[file_name].rootAssembly.instances['Plate_%s'%N_beams].edges.findAt(((Beam_length,Beam_height/float(2),Total_thickness),))
				
				Beam_u3=mdb.models[file_name].rootAssembly.instances['Plate_%s'%N_beams].edges.findAt(((Beam_length,Beam_height/2.,Total_thickness),))
				edges_RBC_u3=Beam_u3
				
				Beam_symm=mdb.models[file_name].rootAssembly.instances['Plate_1'].faces.findAt(((0,Beam_height/2.,Beam_thickness/2.),))
				face_LBC_symm=Beam_symm
				
				Beam_u2=mdb.models[file_name].rootAssembly.instances['Plate_1'].edges.findAt(((Beam_length,0,Beam_thickness/2.),))
				edges_RBC_u2=Beam_u2
				
				if (N_beams>=2):
					for x in range (2,N_beams+1):
						Beam_u2=mdb.models[file_name].rootAssembly.instances['Plate_%s'%x].edges.findAt(((Beam_length,0,Beam_thickness/float(2)+(x-1)*(Beam_thickness+Interlayer_thickness)),))
						edges_RBC_u2=edges_RBC_u2+Beam_u2
					
						Beam_symm=mdb.models[file_name].rootAssembly.instances['Plate_%s'%x].faces.findAt(((0,Beam_height/2.,(x-1)*(Beam_thickness+Interlayer_thickness)+Beam_thickness/float(2)),))
						face_LBC_symm=face_LBC_symm+Beam_symm
				
				
				region_RBC_u3=mdb.models[file_name].rootAssembly.Set(edges=edges_RBC_u3, name='Region_RBC_u3')
				mdb.models[file_name].DisplacementBC(name='BC_RechterSP', createStepName='Initial', region=region_RBC_u3, u3=SET, distributionType=UNIFORM, fieldName='',localCsys=None)
				
				region_RBC_u2=mdb.models[file_name].rootAssembly.Set(edges=edges_RBC_u2, name='Region_RBC_u2')
				mdb.models[file_name].DisplacementBC(name='BC_Rechter_onderrand', createStepName='Initial', region=region_RBC_u2, u2=SET, distributionType=UNIFORM, fieldName='',localCsys=None)
				
				region_LBC=mdb.models[file_name].rootAssembly.Set(faces=face_LBC_symm, name='Region_LBC_symm')
				mdb.models[file_name].XsymmBC(name='BC_Symm', createStepName='Initial', region=region_LBC,localCsys=None)
				
				#_______________Steps and Loads______________#
				
				#Pointload
				Midden_boven=mdb.models[file_name].rootAssembly.instances['Plate_1'].vertices.findAt(((0.,Beam_height,0.),))
				Midden_boven_1=mdb.models[file_name].rootAssembly.instances['Plate_%s'%N_beams].vertices.findAt(((0.,Beam_height,Total_thickness),))
				
				region_load=mdb.models[file_name].rootAssembly.Set(vertices=Midden_boven+Midden_boven_1, name='region_load')
				
				#calculate necessary increments based on buckling load of monolithic beam with comparable properties
				Increments=F_cr(N_beams,Beam_thickness,Beam_height,Beam_length)
				print Increments
				
				mdb.models[file_name].StaticRiksStep(name='GNIA', previous='Initial',maxArcInc=1E3,minArcInc=1E-5,maxNumInc=3000, nlgeom=ON,maxLPF=1.5*int(Increments))
				mdb.models[file_name].ConcentratedForce(name='Concentrated load', createStepName='GNIA', region=region_load, cf2=-0.250, distributionType=UNIFORM, field='', localCsys=None)
				
				SectionRotation_Location=mdb.models[file_name].rootAssembly.instances['Plate_%s'%N_beams].vertices.findAt(((0.,Beam_height,Total_thickness),),((0.,0.,Total_thickness),))
				region_SectionRotation=mdb.models[file_name].rootAssembly.Set(vertices=SectionRotation_Location, name='region_SectionRotation')
			    
				Displacement_Location=mdb.models[file_name].rootAssembly.instances['Plate_%s'%N_beams].vertices.findAt(((0.,Beam_height,Total_thickness),))
				region_Displacement=mdb.models[file_name].rootAssembly.Set(vertices=Displacement_Location, name='region_Displacement')
				Node_number=mdb.models[file_name].rootAssembly.sets['region_Displacement'].nodes[0].label
			    
				#Lateral restraint
				if(lateral_restraint==1):
					Spring_beamend=mdb.models[file_name].parts['Beam'].vertices.findAt(((0,Beam_height,0),),((Beam_length,Beam_height,0),))
					Spring_beamend_subset=mdb.models[file_name].parts['Beam'].Set(vertices=Spring_beamend, name='Spring_beamend_subset')
					Spring_beamend_nodes=mdb.models[file_name].parts['Beam'].sets['Spring_beamend_subset'].nodes
					Spring_beamend_set=mdb.models[file_name].parts['Beam'].Set(nodes=Spring_beamend_nodes, name='Spring_beamend_set')
						
					Spring_edge=mdb.models[file_name].parts['Beam'].edges.findAt(((Beam_length/float(2),Beam_height,0),))
					Spring_subset=mdb.models[file_name].parts['Beam'].Set(edges=Spring_edge, name='Spring_subset')
					Spring_nodes=mdb.models[file_name].parts['Beam'].sets['Spring_subset'].nodes
					Spring_set=mdb.models[file_name].parts['Beam'].Set(nodes=Spring_nodes, name='Spring_set')
					Spring_middle_set=mdb.models[file_name].parts['Beam'].SetByBoolean(name='Spring_middle_set', operation=DIFFERENCE, sets=(Spring_set,Spring_beamend_set,))
				    
					Amount_endnodes=len(Spring_beamend_nodes)
					Amount_middlenodes=len(Spring_nodes)-Amount_endnodes
					Spring_beamend_region=mdb.models[file_name].rootAssembly.instances['Plate_1'].sets['Spring_beamend_set']
					Spring_middle_region=mdb.models[file_name].rootAssembly.instances['Plate_1'].sets['Spring_middle_set']
					mdb.models[file_name].rootAssembly.engineeringFeatures.SpringDashpotToGround(name='Lateral_support_end', region=Spring_beamend_region, orientation=None, dof=3, springBehavior=ON, springStiffness=K_sil/(Amount_middlenodes+1)*0.5, dashpotBehavior=OFF, dashpotCoefficient=0.0)
					mdb.models[file_name].rootAssembly.engineeringFeatures.SpringDashpotToGround(name='Lateral_support_middle', region=Spring_middle_region, orientation=None, dof=3, springBehavior=ON, springStiffness=K_sil/(Amount_middlenodes+1), dashpotBehavior=OFF, dashpotCoefficient=0.0)
		       
				#_______________Job______________#
				del mdb.models[file_name].fieldOutputRequests['F-Output-1']
				del mdb.models[file_name].historyOutputRequests['H-Output-1']
				mdb.models[file_name].FieldOutputRequest(name='GNIA', createStepName='GNIA', variables=('S', 'E',  'U', 'COORD'))
				
				#Add imperfection
				mdb.Job(name='GNA', model=file_name, description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, scratch='D:\Temp', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4, numDomains=4, numGPUs=0)
				#Set temperature
				if(N_beams>=2 and Viscoelastic==1):
					for inst in instances_interlayer:
						region=inst.sets['Interlayerset']
						mdb.models[file_name].Temperature(name='Reference_temperature_%s'%inst.name, createStepName='Initial', region=region,distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(Load_temperature, ))
					for inst in instances_glass:
						region=inst.sets['Beamset']
						mdb.models[file_name].Temperature(name='Reference_temperature_%s'%inst.name, createStepName='Initial', region=region,distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(Load_temperature, ))
		    
				mdb.jobs['GNA'].writeInput()
				
				afterLine='** ----------------------------------------------------------------'
				lines='*Imperfection, input=%s\n'%(impNaam+'.inp')
				addLinesAfter(('GNA.inp'),('GNIA.inp'),lines,afterLine)
				
				t0=time.time()
				path1 = "%s"%dir_temp
				abaqusCall1 = 'abaqus job=GNIA.inp interactive cpus=8'
				runCommand1= 'cmd.exe /c ' + abaqusCall1
				process1 = subprocess.Popen(runCommand1, cwd=path1)
				#~ GNIA_job=mdb.JobFromInputFile(name='GNIA',inputFileName='GNIA', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=FULL, scratch='D:\Temp', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=8, numDomains=8, numGPUs=0)
				file = open("printing.txt", "a")
				file.write("Start infowriting \n")
				file.close()
				#~ GNIA_job.submit()
				file = open("printing.txt", "a")
				file.write("Job has been submitted \n")
				file.close()
				#Determination of exceedance of yield stress
				Buckling_resistance_ANG=0
				Exceeded_ANG=0
				Fail_instance_ANG='no fail'
				Fail_element_ANG='no fail'
				Endframe_ANG=0
				
				Buckling_resistance_HSG=0
				Exceeded_HSG=0
				Fail_instance_HSG='no fail'
				Fail_element_HSG='no fail'
				Endframe_HSG=0
				
				Buckling_resistance_FTG=0
				Exceeded_FTG=0
				Fail_instance_FTG='no fail'
				Fail_element_FTG='no fail'
				Endframe_FTG=0
				
				Stress_at_ip=0
				numFrames=0
				time.sleep(60)
				
				Maximal_stress=0
				file = open("printing.txt", "a")
				file.write("go\n")
				file.close()
				print 'go'
				D_lines=50
				frames=[]
				while (Maximal_stress < Glass_stress('FTG') ):
					file = open("printing.txt", "a")
					fileStat=open('GNIA.sta','r')
					lines=fileStat.readlines()
					n_increments=len(lines)-5
					file.write("Increment_%s \n" %(n_increments))
					if n_increments>=D_lines:
						odb =  ODB_open_safe('GNIA')
						History_regionname=odb.steps['GNIA'].historyRegions.keys(1)[0]
						History_region=odb.steps['GNIA'].historyRegions[History_regionname]
						for instance_name in instances_glass:
							instance_name='%s'%instance_name.name
							instance_name=instance_name.upper()
							stress_elements=odb.rootAssembly.instances['%s'%instance_name].elementSets['STRESS_INCLUSION'].elements
							for x in range (D_lines-50,D_lines):
								frames.append(odb.steps['GNIA'].frames[x])
							for fr in frames:
								for elem in stress_elements:
									Stress=fr.fieldOutputs['S'].getSubset(region=elem).values
									Amount_IP=len(Stress)
									for ip_stress in Stress:
										Stress_at_ip=ip_stress.maxPrincipal
										if(Stress_at_ip>Maximal_stress):
											file = open("printing.txt", "a")
											Maximal_stress=Stress_at_ip
											file.write("Maximal stress %s increment %s \n" %(Maximal_stress,fr.incrementNumber))
											if(Maximal_stress>=Glass_stress('ANG') and Exceeded_ANG<1):
												Buckling_resistance_ANG=History_region.historyOutputs['LPF'].data[fr.incrementNumber][1]
												Endframe_ANG=fr.incrementNumber
												Breakstress_ANG=Maximal_stress
												Fail_instance_ANG=elem.instanceName
												Fail_element_ANG=elem.label
												node_number_ANG=elem.connectivity[0]
												Exceeded_ANG+=1
												print "Exceeded_ANG=" 
												print Exceeded_ANG
											if(Maximal_stress>=Glass_stress('HSG') and Exceeded_HSG<1):
												Buckling_resistance_HSG=History_region.historyOutputs['LPF'].data[fr.incrementNumber][1]
												Endframe_HSG=fr.incrementNumber
												Breakstress_HSG=Maximal_stress
												Fail_instance_HSG=elem.instanceName
												Fail_element_HSG=elem.label
												node_number_HSG=elem.connectivity[0]
												Exceeded_HSG+=1
												print "Exceeded_HSG=" 
												print Exceeded_HSG
											if(Maximal_stress>=Glass_stress('FTG') and Exceeded_FTG<1):
												Buckling_resistance_FTG=History_region.historyOutputs['LPF'].data[fr.incrementNumber][1]
												Endframe_FTG=fr.incrementNumber
												Breakstress_FTG=Maximal_stress
												Fail_instance_FTG=elem.instanceName
												Fail_element_FTG=elem.label
												node_number_FTG=elem.connectivity[0]
												Exceeded_FTG+=1
												print "Exceeded_FTG=" 
												print Exceeded_FTG
												print 'end job'
												file.write("end job \n")
												os.system('taskkill /IM standard.exe /f')
												break
											file.close()
									if(Maximal_stress>=Glass_stress('FTG') ):
										break
								if(Maximal_stress>=Glass_stress('FTG') ):
									break
							if(Maximal_stress>=Glass_stress('FTG') ):
								break
							print Maximal_stress
									
						D_lines+=50
						odb.close()
						frames=[]
					else:
						print 'wachten'
						file.write("wachten\n")
						time.sleep(60)						
				t1=time.time()
				Calctime=t1-t0
				getridkeys=mdb.models[file_name].parts.keys()
				mdb.models[file_name].featureOptions.setValues(autoCaching=OFF)
				for i in range (0,len(getridkeys[i])):
					mdb.models[file_name].parts[getridkeys[i]].clearGeometryCache()
				if (len(mdb.models)>1):
					del mdb.models[file_name]
										  
				##################
				#OUTPUT HANDLING#
				##################
				odb = session.openOdb(name='GNIA.odb')
				session.memoryReductionOptions.setValues(reducedMemoryMode=ON,percentThreshold=50.0)
				session.viewports['Viewport: 1'].setValues(displayedObject=odb)
				session.XYDataFromHistory(name='LPF', odb=odb,outputVariableName='Load proportionality factor: LPF for Whole Model', steps=('GNIA', ), )
				session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', NODAL, ((COMPONENT, 'U3'), )), ), nodeSets=('REGION_DISPLACEMENT', ))
				xy1 = session.xyDataObjects['U:U3 PI: PLATE_%s N: %s'%(N_beams, Node_number)]
				xy2 = session.xyDataObjects['LPF']
				xy3 = combine(xy1, xy2)
				xyp = session.XYPlot('XYPlot-%s'%Plotnr)
				chartName = xyp.charts.keys()[0]
				chart = xyp.charts[chartName]
				c1 = session.Curve(xyData=xy3)
				chart.setValues(curvesToPlot=(c1, ), )
				session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
				session.charts['Chart-%s'%Chartnr].legend.setValues(show=False)
				session.charts['Chart-%s'%Chartnr].gridArea.style.setValues(fill=False)
				session.printOptions.setValues(rendition=GREYSCALE, vpDecorations=OFF)
				session.printToFile(fileName='Load_displacement', format=PNG, canvasObjects=(session.viewports['Viewport: 1'], ))

				session.viewports['Viewport: 1'].setValues(displayedObject=odb)
				session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
				session.viewports['Viewport: 1'].odbDisplay.contourOptions.setValues(maxAutoCompute=OFF, maxValue=Glass_stress('ANG'))
				session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(INVARIANT, 'Max. Principal'), )
				session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(visibleEdges=FEATURE)	    
				session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=Endframe_ANG )
				session.viewports['Viewport: 1'].view.fitView()
				session.printOptions.setValues(rendition=COLOR, vpDecorations=OFF)  
				session.printToFile(fileName=file_name, format=PNG, canvasObjects=(session.viewports['Viewport: 1'], ))
				
				#Output writing
				amount_of_elements=amount_elements(N_beams)
				os.chdir(r"%s"%dir)
				file=open(Analysis_type+'.csv','ab')
				fileResults=csv.writer(file)
				fileResults.writerow([file_name,Calctime,amount_of_elements, Buckling_resistance_ANG,Endframe_ANG,Breakstress_ANG,Fail_instance_ANG,Fail_element_ANG, Estimate_location(node_number_ANG,Fail_instance_ANG,Beam_length,Beam_height),Buckling_resistance_HSG,Endframe_HSG,Breakstress_HSG,Fail_instance_HSG,Fail_element_HSG,Estimate_location(node_number_HSG,Fail_instance_HSG,Beam_length,Beam_height),Buckling_resistance_FTG,Endframe_FTG,Breakstress_FTG,Fail_instance_FTG,Fail_element_FTG,Estimate_location(node_number_FTG,Fail_instance_FTG,Beam_length,Beam_height)])
				getridkeys=session.xyDataObjects.keys() 
				for i in range(0,len(getridkeys)): 
					del session.xyDataObjects[getridkeys[i]] 
				getridkeys=session.xyPlots.keys() 
				for i in range(0,len(getridkeys)): 
					del session.xyPlots[getridkeys[i]] 
				file.close()
				odb.close()
				if (Exceeded_FTG==1):
					mailbody='The %s has finished in %s'%(file_name,Calculation_time(Calctime))
				else:
					mailbody='The %s has not finished correctly'%(file_name)
				send_email('ugent.rekenpc531@gmail.com', 'ir14c531', 'fabrice.kennof@gmail.com', '531:Job finished', 'The job %s has finished in %s'%(file_name,Calculation_time(Calctime)))

	    
	    
	    
#TODO meer algemene benaming (al verbeterd), for lussen voor afmetingen, gelamineerde ligger maken en met input weergeven hoeveel interlayern evt?, efficienter printen,GNIA, problemen met lck.bestand, mesh checken, shear eigenschappen laminaat, mogelijkheden composite layup, BC en Loads inbrengen




