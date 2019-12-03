from math import *
import os
import multiprocessing as mp
from glob import glob

mlfsom_path='/home/peter/Desktop/Summer_Research_2019/mlfsom/'
code_path='/home/peter/Desktop/Summer_Research_2019/code/'

def gaussian(peak=1,sigma=1):
    '''
    Returns a 2d guassian function
    peak - the height of the gaussian
    sigma - the standard deviation of the gaussian
    '''
    return lambda x,y: peak*e**((-x**2-y**2)/(2*sigma**2))

def get_experiment_list(N,start_mos,k_mos,k_cell,k_bfactor,beam_func,frames,frame_rate,radial_symmetry=True,
    angle_slices=1,osc=0.01):
    '''
    Returns a queue of parameters of experiments, weights on the diffraction pattern,
    and frame number, for mlfsom to calculate the spacial and dose dependent diffraction
    patterns for a non-homogenous crystal model given by:
    N by N grid of cells that are centered on the beam
    start_mos - starting value of mosaicity
    k_mos - increase in mosaicity given by mos_0 + k_mos*dose
    k_cell - increase in unit cell dimensions by percentage given by cell_dimenion_0*(1+k_cell*dose)
    k_bfactor - increase in b-factors given by B_0 + k_bfactor*dose
    beam_func - defines spacial distribution of power of beam
    radial_symmetry - if True, assumes beam_func is radially symmetric, if False, does not
    frames - number of frames
    frame_rate - frame rate of experiment
    angle_slices - defines number of angles slices crystal is rotated through in each frames. For a still image set to 1
    osc - angle in degrees that crystal rotates through for each angle slice. For a still image, set osc to 0.01
    Note: dose = frame/frame_rate*beam_func(x,y)
    Returns tuple of list of parameters for experiments to be calculated by mlfsom given by (ID,mos,bfactor_inc,cell_inc,frames,angle_slices,osc)
    and list of weights to each set of experiment on each frame given by (frame,weight,exp_ID)
    '''
    #create list of unique cells and their weights in the crystal
    cell_list=[] #form is (x,y,weight)
    if radial_symmetry:
        #less cells are needed due to symmetry
        if N%2==0:
            for i in range(0,N//2):
                for j in range(0,i+1):
                    if i==j:
                        weight=4
                    else:
                        weight=8
                    cell_list.append((i+.5,j+.5,weight))
        else:
            for i in range(0,N//2+1):
                for j in range(0,i+1):
                    if (i,j)==(0,0):
                        weight=1
                    elif i==j or j==0:
                        weight=4
                    else:
                        weight=8
                    cell_list.append((i,j,weight))
    else:
        #all cells are needed
        if N%2==0:
            for i in range(-N//2,N//2):
                for j in range(-N//2,N//2):
                    cell_list.append((i+.5,j+.5,1))
        else:
            for i in range(-(N//2),N//2+1):
                for j in range(-(N//2),N//2+1):
                    cell_list.append((i,j,1))
    #calculate dose in each cell for each frame and parameters for each cell
    experiment_list=[] #form is (ID,mos,bfactor_inc,cell_inc,frames,osc)
    frame_weights=[] #form is (frame,weight,experiment_ID)
    for frame in range(0,frames):
        exposure=frame/frame_rate
        for cell in cell_list:
            ID=len(experiment_list)
            frame_weights.append((frame+1,cell[2],ID))
            dose=beam_func(cell[0],cell[1])*exposure
            experiment_list.append((ID,start_mos+k_mos*dose,k_bfactor*dose,k_cell*dose,angle_slices,osc))
    #remove repeated experiments and generate list of weights
    i=0
    while i<len(experiment_list):
            experiment=experiment_list[i]
            j=i+1
            while j<len(experiment_list):
                other=experiment_list[j]
                if experiment[1:6]==other[1:6]:
                    experiment_list.pop(j)
                    for index in range(len(frame_weights)):
                        params=frame_weights[index]
                        if params[2]==other[0]:
                            frame_weights[index]=(params[0],params[1],experiment[0])
                else:
                    j+=1
            i+=1
    return (frame_weights,experiment_list)


def get_homogenous_experiment_list(start_mos,k_mos,k_cell,k_bfactor,frames,frame_rate,angle_slices=1,osc=0.01):
    '''
    Returns a queue of parameters of experiments for mlfsom to calculate the dose dependent diffraction
    patterns for a homogenous crystal model given by:
    start_mos - starting value of moasicty
    k_mos - increase in mosaicity given by mos_0 + k_mos*dose
    k_cell - increase in unit cell dimensions by percentage given by cell_dimenion_0*(1+k_cell*dose)
    k_bfactor - increase in b-factors given by B_0 + k_bfactor*dose
    frames - number of frames
    frame_rate - frame rate of experiment (units of 1/time)
    angle_slices - defines number of angles slices crystal is rotated through in each frames. For a still image set to 1
    osc - angle in degrees that crystal rotates through for each angle slice. For a still image, set osc to 0.01
    Note: the dose in each frame is given by frame/frame_rate*1
    Returns list of parameters for experiments to be calculated by mlfsom given by (ID,mos,bfactor_inc,cell_inc,angle_slices,osc)
    '''
    experiment_list=[] #form is (ID,mos,bfactor_inc,cell_inc,frames,osc)
    for frame in range(0,frames):
        exposure=frame/frame_rate
        ID=len(experiment_list)
        dose=1*exposure
        experiment_list.append((ID,start_mos+k_mos*dose,k_bfactor*dose,k_cell*dose,angle_slices,osc))
    return experiment_list

def run_experiments(experiment_list,prefix,threads):
    '''
    Runs each of the experiments and generates the diffraction pattern defined in experiment_list
    Saves images to ~/Desktop/Summer_Research_2019/mlfsom/data/ and necessary temporary files to mlfsom_path/data/tempdata
    experiment_list - list of experiments in form of tuple (ID,mos,bfactor_inc,cell_inc,angle_slices,osc)
    prefix - images are saved as {prefix}###_001.img, where ### is the ID. For instance if prefix is
    mseq and ID is 1, the image will be saved as mseq1_001.img
    threads - number of images that are run at once
    '''
    #copy pdb file to temp file
    prev_cwd=os.getcwd()
    os.chdir(mlfsom_path)
    os.system('cp 1H87.pdb temp.pdb')
    #clear previous temp files
    os.system('rm /tmp/peter/*')
    #set up threading and run experiments
    p=mp.Pool(threads)
    experiments_and_prefix=list(zip(experiment_list,[prefix]*len(experiment_list)))
    for i in range(0,len(experiments_and_prefix),threads):
        #generate input files
        for exp in experiments_and_prefix[i:i+threads]:
            ID=exp[0][0]
            bfactor_inc=exp[0][2]
            cell_inc=exp[0][3]
            print("generating inputs for id="+str(ID))
            os.system('./change_pdb_param.com temp.pdb input'+str(ID)+'.pdb add_cell='+str(cell_inc)+
            ' add_bfactor='+str(bfactor_inc))
            os.system('./ano_sfall.com energy=12660 input'+str(ID)+'.pdb 1.2A solvent_B=35')
            os.system('mv ideal_ano.mtz input'+str(ID)+'.mtz')
        p.map(run_experiment,experiments_and_prefix[i:i+threads])
        #clear temp files and move results
        print("clearing temp files and moving results")
        os.system('mv /tmp/peter/mlfsom*.XYI '+mlfsom_path+'data/tempdata')#change this if you want to save files elsewhere
        os.system('mv /tmp/peter/mlfsom*predin.txt '+mlfsom_path+'data/tempdata')#change this if you want to save files elsewhere
        os.system('rm '+mlfsom_path+'fit2d_*')
        os.system('rm /tmp/peter/*')
        rename_temp_files(experiments_and_prefix[i:i+threads],mlfsom_path+'data/tempdata')
    os.chdir(prev_cwd)

def run_experiment(experiment_and_prefix):
    '''
    Runs mlfsom with parameters defined by experiment and outprefix
    experiment_and_prefix - two tuple, first is experiment parameters given by (ID,mos,bfactor_inc,cell_inc,angle_slices,osc)
    second is prefix
    Images are saved as {prefix}###_001.img, where ### is the ID. For instance if prefix is
    mseq and ID is 1, the image will be saved as mseq1_001.img
    '''
    experiment=experiment_and_prefix[0]
    prefix=experiment_and_prefix[1]
    ID=experiment[0]
    mos=experiment[1]
    bfactor_inc=experiment[2]
    cell_inc=experiment[3]
    frames=experiment[4]
    osc=experiment[5]
    print('running id='+str(ID)+' mos='+str(mos)+' bfactor_inc='+str(bfactor_inc)+' cell_inc='+str(cell_inc)+' osc='+str(osc))
    os.system('./mlfsom.com ~/Desktop/Summer_Research_2019/mlfsom/data/'+prefix+str(ID)+
    '_001.img input'+str(ID)+'.mtz mosaic='+str(mos)+' frames='+str(frames)+' osc='+str(osc)+' id='+str(ID))
    print('----id: '+str(ID)+' is done----')

def spacial_dependent_crystal(N,start_mos,k_mos,k_cell,k_bfactor,beam_func,frames,frame_rate,prefix,angle_slices=1,osc=.01,radial_symmetry=True):
    '''
    Generates files defining the parameters of experiments, weights on the diffraction pattern,
    and frame number, for mlfsom to calculate the spacial and dose dependent diffraction patterns for a non-homogenous crystal model given by:
    N by N grid of cells that are centered on the beam
    start_mos - starting value of mosaicity
    k_mos - increase in mosaicity given by mos_0 + k_mos*dose
    k_cell - increase in unit cell dimensions by percentage given by cell_dimenion_0*(1+k_cell*dose)
    k_bfactor - increase in b-factors given by B_0 + k_bfactor*dose
    beam_func - defines spacial distribution of power of beam, for instance gaussian(peak=1,simga=1)
    radial_symmetry - if True, assumes beam_func is radially symmetric, if False, does not
    frames - number of frames
    frame_rate - frame rate of experiment
    angle_slices - defines number of angles slices crystal is rotated through in each frames. For a still image set to 1
    osc - angle in degrees that crystal rotates through for each angle slice. For a still image, set osc to 0.01
    Files generated are:
    exp-{prefix}-desc.txt - lists the parameters of the experiment
    exp-{prefix}-queue.txt - list of experiments to be run, intended to be where experiments are run from using run_from_exp_queue
    exp-{prefix}-list.txt - back up list of experiments
    Returns tuple of list of parameters for experiments to be calculated by mlfsom given by (ID,mos,bfactor_inc,cell_inc,frames,angle_slices,osc)
    '''
    frame_weights,experiment_list=get_experiment_list(N,start_mos,k_mos,k_cell,k_bfactor,beam_func,frames,frame_rate,radial_symmetry=radial_symmetry,
    angle_slices=angle_slices,osc=osc)
    write_description(N,start_mos,k_mos,k_cell,k_bfactor,beam_func,frames,frame_rate,prefix,angle_slices,osc,frame_weights)
    write_exp_queue_and_list(experiment_list,prefix)
    #run_experiments(experiment_list,prefix,threads)
    return experiment_list


def homogenous_crystal(start_mos,k_mos,k_cell,k_bfactor,frames,frame_rate,prefix,angle_slices=1,osc=.01):
    '''
    Generates files defining the parameters of experiments for mlfsom to calculate the dose dependent diffraction
    patterns for a homogenous crystal model given by:
    start_mos - starting value of moasicty
    k_mos - increase in mosaicity given by mos_0 + k_mos*dose
    k_cell - increase in unit cell dimensions by percentage given by cell_dimenion_0*(1+k_cell*dose)
    k_bfactor - increase in b-factors given by B_0 + k_bfactor*dose
    frames - number of frames
    frame_rate - frame rate of experiment
    angle_slices - defines number of angles slices crystal is rotated through in each frames. For a still image set to 1
    osc - angle in degrees that crystal rotates through for each angle slice. For a still image, set osc to 0.01
    Returns list of parameters for experiments to be calculated by mlfsom given by (ID,mos,bfactor_inc,cell_inc,angle_slices,osc)
    Files generated are:
    exp-{prefix}-desc.txt - lists the parameters of the experiment
    exp-{prefix}-queue.txt - list of experiments to be run, intended to be where experiments are run from using run_from_exp_queue
    exp-{prefix}-list.txt - back up list of experiments
    Returns tuple of list of parameters for experiments to be calculated by mlfsom given by (ID,mos,bfactor_inc,cell_inc,frames,angle_slices,osc)
    '''
    experiment_list=get_homogenous_experiment_list(start_mos,k_mos,k_cell,k_bfactor,frames,frame_rate,angle_slices=angle_slices,osc=osc)
    write_description(1,start_mos,k_mos,k_cell,k_bfactor,"top-hat",frames,frame_rate,prefix,angle_slices,osc,frame_weights=[])
    write_exp_queue_and_list(experiment_list,prefix)
    #run_experiments(experiment_list,prefix,threads)
    return experiment_list

def write_exp_queue_and_list(experiment_list,prefix):
    '''
    Writes list of experiments to exp-{prefix}-queue.txt and exp-{prefix}-list.txt
    experiment_list - list of tuples defing experiment parameters given by (ID,mos,bfactor_inc,cell_inc,angle_slices,osc)
    prefix - name of experiment
    '''
    f=open("exp-"+prefix+"-queue.txt","w")
    flist=open("exp-"+prefix+"-list.txt","w")
    for exp in experiment_list:
        f.write(str(exp)+"\n")
        flist.write(str(exp)+"\n")
    f.close()
    flist.close()

def read_exp_list(file):
    '''
    Reads in experiment list in file
    returns list of tuples defing experiment parameters given by (ID,mos,bfactor_inc,cell_inc,angle_slices,osc)
    '''
    exp_list=[]
    with open(file,"r") as f:
        for i in f:
            exp_list.append(eval(next(f)))
    return exp_list

def run_from_exp_queue(queue_file,N,threads,prefix):
    '''
    Runs and clears top N frames defined in queue_file
    queue_file - file listing frames to be generated by mlfsom
    each frame should be defined by (ID,mos,bfactor_inc,cell_inc,angle_slices,osc)
    Should be the file exp-{prefix}-queue.img
    N - number of frames to be run
    threads - number of frames run at once
    prefix - images saved as {prefix}###_001.img
    Note: frames being run are removed from queue_file, so if they need to be rerun, copy them over from
    exp-{prefix}-list.img
    '''
    exp_list=[]
    with open(queue_file,"r") as f, open("tmp"+queue_file,"w") as out:
        for i in range(N):
            exp_list.append(eval(next(f)))
        for line in f:
            out.write(line)
    os.remove(queue_file)
    os.rename("tmp"+queue_file,queue_file)
    run_experiments(exp_list,prefix,threads)

def write_description(N,start_mos,k_mos,k_cell,k_bfactor,beam_func,frames,frame_rate,prefix,angle_slices,osc,frame_weights):
    f=open("exp-"+prefix+"-desc.txt","w")
    f.write("experiment - "+prefix+"\n")
    f.write("N="+str(N)+"\n")
    f.write("Start Mos="+str(start_mos)+"\n")
    f.write("k_mos="+str(k_mos)+"\n")
    f.write("k_cell="+str(k_cell)+"\n")
    f.write("k_bfactor="+str(k_bfactor)+"\n")
    f.write("beam_func="+str(beam_func)+"\n")
    f.write("frames="+str(frames)+"\n")
    f.write("frame_rate="+str(frame_rate)+"\n")
    f.write("angle_slices="+str(angle_slices)+"\n")
    f.write("osc="+str(osc)+"\n")
    f.write("exp list weights:"+"\n")
    f.write("(frame,weight,experiment_ID)"+"\n")
    for i in frame_weights:
        f.write(str(i)+"\n")
    f.close()

#This does not work correctly
def rename_temp_files(experiments_and_prefix,folder):
    files=glob(folder+"/mlfsom*predin.txt")
    for exp_and_pre in experiments_and_prefix:
        prefix=exp_and_pre[1]
        exp=exp_and_pre[0]
        ID=exp[0]
        mos=exp[1]
        for file in files:
            with open(file) as f:
                lines=f.readlines()
            for l in lines:
                if "ID" in l:
                    fID=l.split()[-1]
                if "MOSAIC" in l:
                    fmos=l.split()[-1]
            if fmos==mos:
                print(ID,mos,fmos,file)

#spacial_dependent_crystal(5,0.05,0.07,0.0153,2.228,gaussian(10/14),14,1,"sp_seq1_")

#To create a queue for a run, run either homogenous_crystal or spacial_dependent_crystal
#Then run, run_from_exp_queue on the queue that you want to run
