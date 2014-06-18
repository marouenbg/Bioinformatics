#! /usr/bin/env python
# -*- coding: utf-8 -*-

#Import Modules
import math
import numpy
import time
import csv
    
def theta_angle(vector1,vector2):
    '''
    Computes theta angle based on the dot product of two vectors
    '''
    cos_theta = numpy.dot(vector1,vector2) / (numpy.linalg.norm(vector1) * numpy.linalg.norm(vector2))
    theta = math.acos(cos_theta)
    return theta
    
def coordinates_two(subunit,starting_residue,pdb):
    '''
    Axis/Gets the coordinates of aminoacids of the first tour of a helix 
    '''
    for i in xrange(len(pdb)):
        if ( pdb[i].startswith("ATOM") and pdb[i][13:16].strip() == "CA" and pdb[i][21] == subunit and int(pdb[i][22:26]) == starting_residue ):
            x1 = float(pdb[i][31:38])
            y1 = float(pdb[i][39:46])
            z1 = float(pdb[i][47:54])
            w1 = pdb[i][17:20]
            vect_1 = numpy.array([x1,y1,z1])
            break
    for j in xrange(i,len(pdb)):               
        if ( pdb[j].startswith("ATOM") and pdb[j][13:16].strip() == "CA" and pdb[j][21] == subunit and int(pdb[j][22:26]) == starting_residue+1 ):
            x2 = float(pdb[j][31:38])
            y2 = float(pdb[j][39:46])
            z2 = float(pdb[j][47:54])
            w2 = pdb[j][17:20]
            vect_2 = numpy.array([x2,y2,z2])
            break
    for k in xrange(j,len(pdb)):
        if ( pdb[k].startswith("ATOM") and pdb[k][13:16].strip() == "CA" and pdb[k][21] == subunit and int(pdb[k][22:26]) == starting_residue+2 ):
            x3 = float(pdb[k][31:38])
            y3 = float(pdb[k][39:46])
            z3 = float(pdb[k][47:54])
            w3 = pdb[k][17:20]
            vect_3 = numpy.array([x3,y3,z3])
            break
    for l in xrange(k,len(pdb)):
        if ( pdb[l].startswith("ATOM") and pdb[l][13:16].strip() == "CA" and pdb[l][21] == subunit and int(pdb[l][22:26]) == starting_residue+3 ):
            x4 = float(pdb[l][31:38])
            y4 = float(pdb[l][39:46])
            z4 = float(pdb[l][47:54])
            w4 = pdb[l][17:20]
            vect_4 = numpy.array([x4,y4,z4])
            break
    aatour=[w1,w2,w3,w4]                  
    return vect_1,vect_2,vect_3,vect_4,aatour
    
def helix_orientation(subunit1,vector1,subunit2,vector2,pdb):
    '''
    Parallel/antiparallel orientation of two helices
    '''
    for i in xrange(len(pdb)):
        if ( pdb[i].startswith("ATOM") and pdb[i][13:16].strip() == "CA" and pdb[i][21] == subunit1 and int(pdb[i][22:26]) == vector1[0] ):
            z1 = float(pdb[i][47:54])            
            break
    for j in xrange(i,len(pdb)):               
        if ( pdb[j].startswith("ATOM") and pdb[j][13:16].strip() == "CA" and pdb[j][21] == subunit1 and int(pdb[j][22:26]) == vector1[1] ):
            z2 = float(pdb[j][47:54])
            break
    for k in xrange(j,len(pdb)):
        if ( pdb[k].startswith("ATOM") and pdb[k][13:16].strip() == "CA" and pdb[k][21] == subunit2 and int(pdb[k][22:26]) == vector2[0] ):
            z3 = float(pdb[k][47:54])
            break
    for l in xrange(k,len(pdb)):
        if ( pdb[l].startswith("ATOM") and pdb[l][13:16].strip() == "CA" and pdb[l][21] == subunit2 and int(pdb[l][22:26]) == vector2[1] ):
            z4 = float(pdb[l][47:54])
            break
    if (z2-z1)*(z4-z3) < 0:
        Or="A"
    else:
        Or="P"
    if (subunit1 == subunit2) :
        Same = "S"
    else:
        Same = "D"                            
    return Or,Same

def coordinates_one(subunit,starting_residue,end_residue,pdb):
    '''
    Midpoint/Gets the coordinates of aminoacids of a helix
    '''
    mid_p=int((end_residue+starting_residue)/2)
    for i in xrange(len(pdb)):
        if ( pdb[i].startswith("ATOM") and pdb[i][13:16].strip() == "CA" and pdb[i][21] == subunit and int(pdb[i][22:26]) == starting_residue ):
            x1 = float(pdb[i][31:38])
            y1 = float(pdb[i][39:46])
            z1 = float(pdb[i][47:54])
            vect_1 = numpy.array([x1,y1,z1])
            break
    for j in xrange(i,len(pdb)):               
        if ( pdb[j].startswith("ATOM") and pdb[j][13:16].strip() == "CA" and pdb[j][21] == subunit and int(pdb[j][22:26]) == mid_p ):
            x2 = float(pdb[j][31:38])
            y2 = float(pdb[j][39:46])
            z2 = float(pdb[j][47:54])
            vect_2 = numpy.array([x2,y2,z2])
            break
    for k in xrange(j,len(pdb)):
        if ( pdb[k].startswith("ATOM") and pdb[k][13:16].strip() == "CA" and pdb[k][21] == subunit and int(pdb[k][22:26]) == end_residue ):
            x3 = float(pdb[k][31:38])
            y3 = float(pdb[k][39:46])
            z3 = float(pdb[k][47:54])
            vect_3 = numpy.array([x3,y3,z3])
            break
    return vect_1,vect_2,vect_3

def helix_length_angs(subunit,helix_start,helix_end,pdb):
    '''
    Computes helix length in Angstroms and the helix center without 
    taking the kinks into account
    '''
    vect_coor= coordinates_one(subunit,helix_start,helix_end,pdb)
    start_vector = vect_coor[0]
    end_vector = vect_coor[2]
    helix_mid= vect_coor[1]
    vect_coor=[]
    helix_length_angs = math.sqrt( (end_vector[0] - start_vector[0])**2 + (end_vector[1] - start_vector[1])**2 + (end_vector[2] - start_vector[2])**2 )
    helix_center = [(end_vector[0] + start_vector[0])/2,(end_vector[1] + start_vector[1])/2,(end_vector[2] + start_vector[2])/2]
    return helix_length_angs,helix_center,helix_mid

def axis(subunit,starting_residue,pdb):
    '''
    Does the first step of axis calculation, computes local vector and AA seuence of a tour
    '''	
    t,t1=[0,0,0],[0,0,0]
    result = coordinates_two(subunit,starting_residue,pdb)
    a = result[0]#Coordinates of each AA of the tour
    b = result[1]
    c = result[2]
    d = result[3]
    e = result[4]#The sequence of the tour
    result=[]
    vec1 = [i - j for i, j in zip(b, a)]
    vec2 = [i - j for i, j in zip(c, b)]
    vec3 = [i - j for i, j in zip(d, c)]
    vec4 = [i - j for i, j in zip(vec1, vec2)]
    vec5 = [i - j for i, j in zip(vec2, vec3)] 
    t[0] = (vec4[1] * vec5[2]) - (vec4[2] * vec5[1])
    t[1] = (vec4[2] * vec5[0]) - (vec4[0] * vec5[2])
    t[2] = (vec4[0] * vec5[1]) - (vec4[1] * vec5[0])
    mag=numpy.linalg.norm(t)
    t1[0]=t[0]/mag
    t1[1]=t[1]/mag
    t1[2]=t[2]/mag
    return t1,e
    
def axis_final(unit,h_start,h_end,pdb):
    '''
    Computes final axis, mean of all local vectors and final AA sequence
    '''	
    vect_t=[0,0,0]
    count=0
    angle_list,first_tour,last_tour=[],[],[]
    for i in xrange(h_start,h_end-3):
        aa=axis(unit,i,pdb)
        a=aa[0]
        j=i-h_start
        vect_t[0]=a[0]+vect_t[0]
        vect_t[1]=a[1]+vect_t[1]
        vect_t[2]=a[2]+vect_t[2]
        if i == h_start :
            first_tour.append(aa[1])
        if i == h_end-4 :
            last_tour.append(axis(unit,i+1,pdb)[1])        
        if j==0:
            b=a
        if (j%3==0) and (j != 0):
            angle = a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
            b=a
            anglef=math.degrees(math.acos(angle))
            if anglef>=20:
                count=count+1    
            angle_list.append(anglef)
    vect_t[0]=vect_t[0]/(h_end-h_start-2)
    vect_t[1]=vect_t[1]/(h_end-h_start-2)
    vect_t[2]=vect_t[2]/(h_end-h_start-2)
    return vect_t,count,max(angle_list),first_tour,last_tour

def tilt_angle(helx_ax):
    '''
    Computes tilt angle in regards to the bilayer normal
    '''
    prot_axis=[0,0,1]
    angle = math.degrees(theta_angle(helx_ax,prot_axis))
    if angle>90:
        angle= 180-angle
    return angle
    
def helices_ininteraction(subunit_vector,subunit_helix,pdb_coordinates):
    '''
    Helices are considered in interaction if three or more residues
    were in interaction(chothia et al, 1981), residues in contact
    are defined by a threshold of 10A betwenn their respective CA
    '''
    scale = {"ALA":0.31,"CYS":1.54,"ASP":-0.77,"GLU":-0.64,"PHE":1.79,\
    "GLY":0,"HIS":0.13,"ILE":1.8,"LYS":-0.99,"LEU":1.7,"MET":1.23,\
    "ASN":-0.6,"PRO":0.72,"GLN":-0.22,"ARG":-1.01,"SER":-0.04,"THR":0.26,\
    "VAL":1.22,"TRP":2.25,"TYR":0.96}
    CA_matrix,moment_list,m_hydro_list = [],[],[]
    for i in xrange(len(subunit_vector)):
        CA_subunit,m_hydro,hmoment = [],[],[]
        for j in xrange(len(subunit_helix[i])):
            CA_array,mean_hydro,f_term,s_term = [],0,0,0
            for residue_id in xrange(subunit_helix[i][j][0],subunit_helix[i][j][1]+1):
                for line in pdb_coordinates:        
                    if ( line.startswith("ATOM") and line[13:16].strip() == "CA" and line[21] == subunit_vector[i] and int(line[22:26]) == residue_id ):
                        x1 = float(line[31:38])
                        y1 = float(line[39:46])
                        z1 = float(line[47:54])
                        CA1 = [x1,y1,z1]
                        CA_array.append(CA1)
                        CA1 = []
                        H = scale[line[17:20]]
                        L = subunit_helix[i][j][1]-subunit_helix[i][j][0]+1
                        #mean hydrophobicity
                        mean_hydro = mean_hydro+(H/L)
                        #mean amphiphatic moment
                        n=residue_id-subunit_helix[i][j][0]
                        f_term = f_term+H*math.sin(n*1.7453)
                        s_term = s_term+H*math.cos(n*1.7453)
            hydro_moment=math.sqrt(f_term**2+s_term**2)/L
            hmoment.append(hydro_moment)
            hydro_moment=0
            m_hydro.append(mean_hydro)
            mean_hydro=0
            CA_subunit.append(CA_array)
            CA_array=[]            
        CA_matrix.append(CA_subunit)
        CA_subunit=[]
        moment_list.append(hmoment)
        hmoment=[]
        m_hydro_list.append(m_hydro)
        m_hydro=[]
    int_list,int_list2=[],[]
    closest_app,c_vector,h_orient,Hel_S=[],[],[],[]
    for i in xrange(len(CA_matrix)):
        for j in xrange(len(CA_matrix[i])):
            for k in xrange(len(CA_matrix[i][j])):
                if [i+1,j+1,k+1] == [len(CA_matrix),len(CA_matrix[i]),len(CA_matrix[i][j])]:
                    break
                A = numpy.array((CA_matrix[i][j][k]))
                for l in xrange(i,len(CA_matrix)):
                    if i == l:
                        ss = j+1
                    else:
                        ss = 0        
                    for m in xrange(ss,len(CA_matrix[l])):
                        for n in xrange(len(CA_matrix[l][m])): 
                            B = numpy.array((CA_matrix[l][m][n]))
                            aa12=helix_length_angs(subunit_vector[i],subunit_helix[i][j][0],subunit_helix[i][j][1],pdb_coordinates)[2]
                            bb12=helix_length_angs(subunit_vector[l],subunit_helix[l][m][0],subunit_helix[l][m][1],pdb_coordinates)[2]
                            if (numpy.linalg.norm(A-B) <= 10) and (numpy.linalg.norm(aa12-bb12) <= 25) :
                                int_list.append([subunit_vector[i],j,subunit_vector[l],m])
                                int_list2.append([subunit_vector[i]+str(j),subunit_vector[l]+str(m)])
                                closest_app.append(numpy.linalg.norm(A-B))
                                v = [o - r for o, r in zip(A, B)]
                                c_vector.append(v)
                                v=0
				dresult = helix_orientation(subunit_vector[i],subunit_helix[i][j],subunit_vector[l],subunit_helix[l][m],pdb_coordinates)
                                h_orient.append(dresult[0])
                                Hel_S.append(dresult[1])
                                continue
    CA_matrix,h_or_final=[],[]
    int_helices,h_or,int_helices2=[],[],[]
    Hel_SS,closest,c_vect_final=[],[],[]
    final_count,num_inter,cvect=[],[],[]                            
    for i in xrange(len(int_list)):
        count = 0
        c = [closest_app[i]]
        cvect = [c_vector[i]]
        for j in xrange(i+1,len(int_list)):
            if int_list[j]==int_list[i]:
                count = count+1
                c.append(closest_app[j])
                cvect.append(c_vector[j])
        if count>=3:
            final_count.append(count)
            int_helices.append(int_list[i])
            int_helices2.append(int_list2[i])
            closest.append(min(c))
            ind = numpy.argmin(c)
            c_vect_final.append(cvect[ind])
            h_or.append(h_orient[i])
            Hel_SS.append(Hel_S[i])    
    interacting_helices,interacting_helices2=[],[]
    c_approach,cvectf,cvectfnorm,Hel_SSD = [],[],[],[]
    for i in xrange(len(int_helices)):
        if int_helices[i] not in interacting_helices:
            interacting_helices.append(int_helices[i])
            interacting_helices2.append(int_helices2[i])
            c_approach.append(closest[i])
            num_inter.append(final_count[i])
            cvectf.append(c_vect_final[i])
            h_or_final.append(h_or[i])
            Hel_SSD.append(Hel_SS[i])
            cvectfnorm.append(numpy.linalg.norm(c_vect_final[i]))
    return interacting_helices,c_approach,num_inter,moment_list,m_hydro_list,cvectf,h_or_final,interacting_helices2,cvectfnorm,Hel_SSD
    
def tilt_length_stats(subunit_vector,subunit_helix,pdb_coordinates):
    '''
    Gives out all statistics related to tilt angles, length of helices, number of kinks and maximum kink angles
    '''
    helix_length_tot=[]
    helix_center_list = []
    for k in xrange(len(subunit_vector)):
        centers = []
        for i in subunit_helix[k]:
            helix_start = i[0]
            helix_end = i[1]
            helix_length_tot.append(helix_length_angs(subunit_vector[k],helix_start,helix_end,pdb_coordinates)[0])
            centers.append(helix_length_angs(subunit_vector[k],helix_start,helix_end,pdb_coordinates)[1])
        helix_center_list.append(centers)
    moylen=numpy.mean(helix_length_tot)
    stdlen=numpy.std(helix_length_tot)
    helix_tilt_tot,kangles=[],[]
    first_tours,last_tours,kinks_f=[],[],[]
    for k in xrange(len(subunit_vector)):
        kinks=[]        
        for i in xrange(len(subunit_helix[k])):
            first_tour,last_tour=[],[]
            helix_start = subunit_helix[k][i][0]
            helix_end = subunit_helix[k][i][1]
            z=axis_final(subunit_vector[k],helix_start,helix_end,pdb_coordinates)
            helix_tilt_tot.append(tilt_angle(z[0]))
            kinks.append(z[1])
            kangles.append(z[2])
            first_tours.append(sum(z[3],[]))
            last_tours.append(sum(z[4],[]))
        kinks_f.append(kinks)    
    dist_aa_tot=[]
    for i in subunit_helix:
        for j in i:
            dist_aa=j[1]-j[0]+1
            dist_aa_tot.append(dist_aa)
    distmoyaa,diststdaa=numpy.mean(dist_aa_tot),numpy.std(dist_aa_tot)
    tiltlen,std_len = numpy.mean(helix_tilt_tot),numpy.std(helix_tilt_tot)
    return helix_tilt_tot,helix_length_tot,dist_aa_tot,moylen,stdlen,tiltlen,std_len,distmoyaa,diststdaa,helix_center_list,kangles,kinks_f,first_tours,last_tours

def packing_stats(subunit_vector,subunit_helix,pdb_coordinates):
    '''
    Gives all statistics regarding packing angles
    '''
    result=helices_ininteraction(subunit_vector,subunit_helix,pdb_coordinates)
    h_or=result[6]
    list_interhel=result[0]
    clos_vect=result[5]
    result=[]
    packing_angles,Or = [],[]
    for i in xrange(len(list_interhel)):
        for j in xrange(len(subunit_vector)):
            if subunit_vector[j]==list_interhel[i][0]:
                hel1 = axis_final(list_interhel[i][0],subunit_helix[j][list_interhel[i][1]][0],subunit_helix[j][list_interhel[i][1]][1],pdb_coordinates)[0]
            if subunit_vector[j]==list_interhel[i][2]:
                hel2 = axis_final(list_interhel[i][2],subunit_helix[j][list_interhel[i][3]][0],subunit_helix[j][list_interhel[i][3]][1],pdb_coordinates)[0]
        n1 = numpy.cross(clos_vect[i],hel1)
        n2 = numpy.cross(clos_vect[i],hel2)
        ref = numpy.dot(n1,n2)
        omega = math.degrees(numpy.arccos(ref / (numpy.linalg.norm(n1) * numpy.linalg.norm(n2))))
        volume=numpy.dot(n2,hel1);
        if volume<0:
            omega=-omega
        if omega<-90:
            omega = 180 + omega
        if omega>90:
            omega = -(180 - omega)
        packing_angles.append(omega)
    clos_vect=[]  
    return packing_angles,list_interhel,h_or


def neighbor_hel(subunit_vector,subunit_helix,pdb_coordinates):
    '''
    Computes the percentage of sequence neighbor helices that are in intereaction
    '''
    h_seq=[]
    l=len(sum(subunit_helix,[]))
    n_int=0
    hel=helices_ininteraction(subunit_vector,subunit_helix,pdb_coordinates)[0]
    for i in hel:
        for j in xrange(len(subunit_vector)):
            if i[0]==subunit_vector[j]:
                a=len(subunit_helix[j])
            if i[2]==subunit_vector[j]:
                b=len(subunit_helix[j])        
        if i[0]==i[2]:
            if (i[3]-i[1]==-1) or (i[3]-i[1]==1):
                n_int = n_int+1
                h_seq.append("N")
            else:
                h_seq.append("NOTN")
        else:
            if (i[1]==0 and i[3]==b-1) or (i[3]==0 and i[1]==a-1):
                n_int = n_int+1
                h_seq.append("N")
            else:
                h_seq.append("NOTN")
    neighbor_helix_int=float(n_int)/float(l-1)
    return neighbor_helix_int,h_seq
    
def connexion_loop(subunit_vector,subunit_helix,pdb_coordinates):
    '''
    Computes the length of the connexion loop between interacting helices
    '''
    hel=helices_ininteraction(subunit_vector,subunit_helix,pdb_coordinates)[0]
    cnx_list=[]
    for i in hel:
        for j in xrange(len(subunit_vector)):
            if i[0]==subunit_vector[j]:
                aa=subunit_helix[j][i[1]][0]
                bb=subunit_helix[j][i[1]][1]
                for k in xrange(len(subunit_vector)):
                    if i[2]==subunit_vector[k]:
                        cc=subunit_helix[k][i[3]][0]
                        dd=subunit_helix[k][i[3]][1]
        cnx_length=cc-bb        
        if i[0]==i[2] :        
            if cc-bb>0:
                cnx_list.append(cc-bb)
            else:
                cnx_list.append(aa-dd)
        if i[0]!=i[2] :
            for m in xrange(len(pdb_coordinates)):
                if ( pdb_coordinates[m].startswith("ATOM") and pdb_coordinates[m][21] == i[0] and int(pdb_coordinates[m][22:26]) == bb ):
                    limit1=m
                    break
            for zz in xrange(len(pdb_coordinates)):
                if ( pdb_coordinates[zz].startswith("ATOM") and pdb_coordinates[zz][21] == i[2] and int(pdb_coordinates[zz][22:26]) == cc ):
                    limit2=zz
                    break
            for l in xrange(min(limit1,limit2),max(limit1,limit2)):
		if (pdb_coordinates[l].startswith("HETATM")):
		    pdb_coordinates[l] = pdb_coordinates[l-1]
                if (pdb_coordinates[l][21] != pdb_coordinates[l-1][21]) :            
                    cnx_length=cnx_length+int(pdb_coordinates[l-1][22:26])            
            cnx_list.append(cnx_length)
    return cnx_list
    
def neighbor_helices(subunit_vector,subunit_helix,pdb_coordinates):
    '''
    Computes the number of helices that each helix interacts with
    '''
    a=sum(helices_ininteraction(subunit_vector,subunit_helix,pdb_coordinates)[7],[])
    freq_list,freq_list2=[],[]
    for i in xrange(len(subunit_vector)):
        for j in xrange(len(subunit_helix[i])):
            count=0
            for k in a:
                if k==subunit_vector[i]+str(j):
                    count=count+1
            freq_list.append(count)
        freq_list2.append(freq_list)            
    return freq_list2

def pack_st(subunit_vector,subunit_helix,pdb_coordinates):
    '''
    Eisenberg Amphiphatic mean moment of hydrophobicity
    '''    
    a=sum(helices_ininteraction(subunit_vector,subunit_helix,pdb_coordinates)[7],[])
    c=helices_ininteraction(subunit_vector,subunit_helix,pdb_coordinates)[4]
    kink_list=tilt_length_stats(subunit_vector,subunit_helix,pdb_coordinates)[11]    
    firstkink,secondkink,firsthydro,secondhydro=[],[],[],[]    
    for k in xrange(len(a)):
        for i in xrange(len(subunit_vector)):
            for j in xrange(len(subunit_helix[i])):
                if a[k]==subunit_vector[i]+str(j):
                    if k%2==0:                    
                        firstkink.append(kink_list[i][j])
                        firsthydro.append(c[i][j])
                    if k%2==1:
                        secondkink.append(kink_list[i][j])
                        secondhydro.append(c[i][j])    
    return firstkink,secondkink,firsthydro,secondhydro

def main():
    start=time.time()
    
    configfile = open("config")
    prot_list = configfile.readline()[:-1].split(",")
    subunit_l = configfile.readline()[:-1]
    helix_l = configfile.readline()[:-1]
    subunit_list,helix_list,helix_l1,helix_l2 = [],[],[],[]
    for i in xrange(len(subunit_l.split(";"))):
        subunit_list.append(subunit_l.split(";")[i].split(","))
    for i in xrange(len(helix_l.split(";"))):
        helix_l2=[]
        for j in xrange(len(helix_l.split(";")[i].split("/"))):
            helix_l1=[]
            for k in xrange(len(helix_l.split(";")[i].split("/")[j].split(":"))):
                helix_l1.append(helix_l.split(";")[i].split("/")[j].split(":")[k].split(","))
                for l in xrange(len(helix_l1)):
                    helix_l1[l][0]=int(helix_l1[l][0])
                    helix_l1[l][1]=int(helix_l1[l][1])
            helix_l2.append(helix_l1)
        helix_list.append(helix_l2) 
   
    nh,hmoy,ap,contacts,closest,loop,first,second,fhydro,shydro,kink_angle,n_kinks,tilt=[],[],[],[],[],[],[],[],[],[],[],[],[]    
    SD,nei_hel,heisenberg,packing,lengtha,lengthaa,first_tour,last_tour,neigh_seq=[],[],[],[],[],[],[],[],[]   
    for i in xrange(len(prot_list)):
        fichier = open("%s.pdb"%prot_list[i])
        subunit_vector = subunit_list[i]
        subunit_helix = helix_list[i]

        pdb_coordinates = fichier.readlines()
        fichier.close()

    #print helices_ininteraction(subunit_vector,subunit_helix,pdb_coordinates)[2]
        #Closest approach
        a=packing_stats(subunit_vector,subunit_helix,pdb_coordinates)
        #packing angles           
        b=tilt_length_stats(subunit_vector,subunit_helix,pdb_coordinates)
        c=helices_ininteraction(subunit_vector,subunit_helix,pdb_coordinates)
        d=neighbor_helices(subunit_vector,subunit_helix,pdb_coordinates)
        e=pack_st(subunit_vector,subunit_helix,pdb_coordinates)
        f=neighbor_hel(subunit_vector,subunit_helix,pdb_coordinates)
    #Outpout file
        neigh_seq.append(f[1])
        kink_angle.append(sum(b[11],[]))
        tilt.append(b[0])
        lengtha.append(b[1])
        lengthaa.append(b[2])
        last_tour.append(b[13])
        first_tour.append(b[12])
        nh.append(d[0])
        n_kinks.append(b[10])
        packing.append(a[0])
        ap.append(a[2])
        contacts.append(c[2])
        closest.append(c[8])
        SD.append(c[9])
        loop.append(connexion_loop(subunit_vector,subunit_helix,pdb_coordinates))
        first.append(e[0])
        second.append(e[1])
        fhydro.append(e[2])
        shydro.append(e[3])
        hmoy.append(sum(c[3], []))
        heisenberg.append(sum(c[4], []))
        nei_hel.append(f[0])
        print("%s done"%prot_list[i])
    myfile = open("Packing_angles.csv",'wb')#Packing angles stats
    wr = csv.writer(myfile,delimiter=',' ,quoting=csv.QUOTE_ALL)
    wr.writerow(sum(packing, [])) #Packing angles
    wr.writerow(sum(ap, [])) #A/P orientation
    wr.writerow(sum(contacts, [])) #Number of contacts between helices
    wr.writerow(sum(closest, [])) #Distance of closest approach
    wr.writerow(sum(loop, [])) #Length of connexion loop     
    wr.writerow(sum(first, [])) #Number of kinks in the first helix
    wr.writerow(sum(second, [])) #Number of kinks in the second helix
    wr.writerow(sum(fhydro, [])) #Eisenberg hydrophobicity in the first helix
    wr.writerow(sum(shydro, [])) #Eisenberg hydrophobicity in the second helix
    wr.writerow([numpy.mean(nei_hel)]) #Percentage of interacting helices that are sequence neighbors        	
    wr.writerow(sum(SD, []))#Interchain interaction
    wr.writerow(sum(neigh_seq, []))#Helices neighbors in sequence
    myfile.close()
    myfile2 = open("Tilt_angles.csv",'wb')#Tilt angles stats
    wr = csv.writer(myfile2,delimiter=',' ,quoting=csv.QUOTE_ALL)
    wr.writerow(sum(tilt, [])) #Tilt angles
    wr.writerow([b[5]]) #Mean
    wr.writerow([b[6]]) #SD
    wr.writerow(sum(lengtha, [])) #Length in AA residus
    wr.writerow([b[3]]) #Mean
    wr.writerow([b[4]]) #SD
    wr.writerow(sum(lengthaa, [])) #Length Angstrom
    wr.writerow([b[7]]) #Mean
    wr.writerow([b[8]]) #SD
    wr.writerow(sum(hmoy, []))#Mean hydrophobicity
    wr.writerow(sum(heisenberg, []))#Eisenberg hydrophobicity
    wr.writerow(sum(nh, [])) #Neighbor helices
    wr.writerow(sum(n_kinks, [])) #largest kink angle
    wr.writerow(sum(kink_angle, [])) #Number of kinks per helix
    wr.writerow(sum(first_tour, []))#Sequence of first tour
    wr.writerow(sum(last_tour, []))#Sequence of last tour
    myfile2.close()
    print time.time()-start
if __name__ == '__main__':
    main()
