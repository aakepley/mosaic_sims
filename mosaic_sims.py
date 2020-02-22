import os

def runTests(vis, casa1, casa2, image_params={'robust':'0.5',
                                              'imsize':'[140,140]',
                                              'cell':'2.4arcsec',
                                              'phasecenter':'ICRS 05:00:00.0 -053.00.00.0'}):

    '''
    setup the appropriate tests
    '''
   
    print('Estimating the RMS')
    estimateRMS(vis,casa1, image_params)
    os.system(casa1['path']+'/bin/casa -c estimate.py')

    print('Creating a mask')
    createMask(vis,casa1, image_params)
    os.system(casa1['path'] + '/bin/casa -c mask.py')

    print('Cleaning the first image')
    createCleanImage(vis,casa1,image_params)
    os.system(casa1['path']+'/bin/casa' + ' -c tclean_'+casa1['name']+'.py')

    print ('Cleaning the second image')
    createCleanImage(vis,casa2, image_params)
    os.system(casa2['path']+'/bin/casa' + ' -c tclean_'+casa2['name']+'.py')

    print('Doing Analysis')
    createAnalysis(vis,casa1,casa2)
    os.system(casa1['path']+'/bin/casa' + ' -c analysis.py')


def estimateRMS(vis,casa,image_params):

    '''
    create a dirty image and use that to estimate the RMS noise
    '''


    outStr = '''
import os 

vis = '{0:s}'
robust = {1:s}
imsize = {2:s}
cell = '{3:s}'
phasecenter = '{4:s}'
    
imagename=os.path.basename(vis).replace('noise.aca.cycle6.','')+'.dirty_'+'{5:s}'
os.system('rm -rf '+imagename+'.*')
tclean(vis=vis,
       imsize=imsize,cell=cell,
       imagename=imagename,
       specmode='mfs',
       phasecenter=phasecenter,mosweight=True,
       gridder='mosaic', pblimit=0.2, deconvolver='hogbom',
       restoration=True,pbcor=True,
       weighting='briggs', robust=robust, niter=0, threshold='0.0Jy',
       interactive=0, savemodel='none', parallel=False)

dirtyrms = 0.5 * 1.4825 * aU.imageImstat(img=imagename+'.image',
                                         keyword='medabsdevmed',pb=imagename+'.pb',
                                         pbInner = 0.35, pbOuter = 0.2)

fout = open('rms.txt','w')
fout.write(str(dirtyrms))
fout.close()
    '''.format(vis,
               image_params['robust'],
               image_params['imsize'],
               image_params['cell'],
               image_params['phasecenter'],
               casa['name'])

    print(outStr)

    fout = open('estimate.py','w')
    fout.write("##Run in "+casa['path']+"/bin/casa \n")
    fout.write(outStr)
    fout.close()

def createMask(vis,casa,image_params):

    '''
    use auto-multithresh to create a mask for the image
    '''


    f = open('rms.txt','r')
    dirtyrms = f.readline()
    f.close()

    if 'sidelobethreshold' in image_params:
        sidelobethreshold = image_params['sidelobethreshold']
    else:
        sidelobethreshold = '1.25'

    if 'cyclefactor' in image_params:
        cyclefactor = image_params['cyclefactor']
    else:
        cyclefactor='1.0'

    outStr = '''
vis = '{0:s}'
robust = {1:s}
imsize = {2:s}
cell = '{3:s}'
phasecenter = '{4:s}'
sidelobethreshold = {7:s}

imagename=os.path.basename(vis).replace('noise.aca.cycle6.','')+'.clean_'+'{5:s}'
os.system('rm -rf '+imagename+'.*')
tclean(vis=vis,
       imsize=imsize,cell=cell,
       imagename=imagename,
       specmode='mfs',
       phasecenter=phasecenter,mosweight=True,
       gridder='mosaic', pblimit=0.2, deconvolver='hogbom',
       restoration=True,pbcor=True,
       weighting='briggs', robust=robust, niter=100000, threshold='{6:s}Jy',
       interactive=0, savemodel='none', parallel=False,
       usemask='auto-multithresh', sidelobethreshold=sidelobethreshold,smoothfactor=0.7,
       noisethreshold=5.0, lownoisethreshold=4.00, negativethreshold=0.0,
       minbeamfrac=0.3, growiterations=75, dogrowprune=True,
       minpercentchange=1.0,cyclefactor={8:s})

os.system('cp -r '+imagename+'.mask final.mask') 
    '''.format(vis, #0
               image_params['robust'], #1
               image_params['imsize'], #2
               image_params['cell'], #3
               image_params['phasecenter'], #4
               casa['name'], #5
               dirtyrms,#6
               sidelobethreshold, #7
               cyclefactor) #8

    print(outStr)

    fout = open('mask.py','w')
    fout.write("##Run in "+casa['path']+"/bin/casa \n")
    fout.write(outStr)
    fout.close()

def createCleanImage(vis,casa,image_params):

    '''
    clean the image using the previously determined mask and threshold.
    '''

    f = open('rms.txt','r')
    dirtyrms = f.readline()
    f.close()

    if 'cyclefactor' in image_params:
        cyclefactor = image_params['cyclefactor']
    else:
        cyclefactor='1.0'

    outStr = '''
vis = '{0:s}'
robust = {1:s}
imsize = {2:s}
cell = '{3:s}'
phasecenter = '{4:s}'

imagename=os.path.basename(vis).replace('noise.aca.cycle6.','')+'.dirty_'+'{5:s}'
os.system('rm -rf '+imagename+'.*')
tclean(vis=vis,
       imsize=imsize,cell=cell,
       imagename=imagename,
       specmode='mfs',
       phasecenter=phasecenter,mosweight=True,
       gridder='mosaic', pblimit=0.2, deconvolver='hogbom',
       restoration=True,pbcor=True,
       weighting='briggs', robust=robust, niter=0, threshold='0.0Jy',
       interactive=0, savemodel='none', parallel=False)


imagename=os.path.basename(vis).replace('noise.aca.cycle6.','')+'.clean_'+'{5:s}'
os.system('rm -rf '+imagename+'.*')
tclean(vis=vis,
       imsize=imsize,cell=cell,
       imagename=imagename,
       specmode='mfs',
       phasecenter=phasecenter,mosweight=True,
       gridder='mosaic', pblimit=0.2, deconvolver='hogbom',
       restoration=True,pbcor=True,
       weighting='briggs', robust=robust, niter=100000, threshold='{6:s}Jy',
       interactive=0, savemodel='none', parallel=False,
       mask='final.mask',cyclefactor={7:s})
    
    '''.format(vis, #0
               image_params['robust'], #1
               image_params['imsize'], #2
               image_params['cell'], #3
               image_params['phasecenter'], #4
               casa['name'],dirtyrms, #5, #6
               cyclefactor) #7

    filename = 'tclean_'+casa['name']+'.py'

    fout = open(filename,'w')
    fout.write("##Run in "+casa['path']+"/bin/casa \n")
    fout.write(outStr)
    fout.close()




def createAnalysis(vis,casa1,casa2):

    '''
    compare the images produced by two different versions of casa
    '''

    import re

    tmp = re.search('_(\d\d)x(\d\d)',vis)
    simsize = [int(tmp[1]),int(tmp[2])]

    array = int(re.search('_(\d\d)m_',vis)[1])
    

    outStr = '''
import analyzemsimage as ami
import csv

vis = '{0:s}'
img0=os.path.basename(vis).replace('noise.aca.cycle6.','')+'.clean_'+'{1:s}'+'.image.pbcor'
img1=os.path.basename(vis).replace('noise.aca.cycle6.','')+'.clean_'+'{2:s}'+'.image.pbcor'
    
rms0 = ami.runImstat(Image=img0,
                     PB=img0.replace('image.pbcor','pb'),
                     Mask=img0.replace('image.pbcor','mask'),
                     innerAnnulusLevel=0.3, level=0.2)['rms'][0]

rms1 = ami.runImstat(Image=img1,
                     PB=img0.replace('image.pbcor','pb'),
                     Mask=img1.replace('image.pbcor','mask'),
                     innerAnnulusLevel=0.3, level=0.2)['rms'][0]

# Do S/N cut of pdiff based on 5.6 rms noise 
cut = 7 * rms1
shortvis = vis.replace('noise.aca.cycle6.','')
pdiff = img0.replace('image.pbcor','')+'vs.'+img1.replace(shortvis+'.clean_','')+'_pdiff'

os.system('rm -rf '+pdiff)
immath(imagename=[img0,img1],
       expr='100*(IM0-IM1)/IM0',
       mask=img0+'>'+str(cut)+' && '+img1+'>'+str(cut),
       outfile=pdiff)

diff = img0.replace('image.pbcor','')+'vs.'+img1.replace(shortvis+'.clean_','')+'_diff'
os.system('rm -rf '+diff)
immath(imagename=[img0,img1],
        expr='(IM0-IM1)',
        outfile=diff)

statspdiff = ami.runImstat(Image=pdiff,
                           PB=img0.replace('image.pbcor','pb'),
                           innerAnnulusLevel=1.0, level=0.5)
outstatspdiff = ami.runImstat(Image=pdiff,
                              PB=img0.replace('image.pbcor','pb'),
                              innerAnnulusLevel=0.5, level=0.2)
statsdiff = ami.runImstat(Image=diff,
                          PB=img0.replace('image.pbcor','pb'),
                          innerAnnulusLevel=1.0, level=0.5)
outstatsdiff = ami.runImstat(Image=diff,
                             PB=img0.replace('image.pbcor','pb'),
                             innerAnnulusLevel=0.5, level=0.2)

aU.imviewFields(imgs=[img0,img1],
                minIntensity=-0.001, maxIntensity=0.21,scaling=-0.75,
                contourImage=img0.replace('image.pbcor','pb'), 
                levels=[0.5], unit=1.0, contourThickness=2) 

aU.imviewFields(imgs=[diff],
                scaling=-0.75,
                contourImage=img0.replace('image.pbcor','pb'), 
                levels=[0.5], unit=1.0, contourThickness=2)           

# remove Jy/beam from pdiff so it won't show on figure, can't use non-casa unit.
imhead(imagename=pdiff,mode='put',hdkey='bunit',hdvalue='')

aU.imviewFields(imgs=[pdiff],
                contourImage=img0.replace('image.pbcor','pb'), 
                levels=[0.5], unit=1.0, contourThickness=2)

# Analyze difference in .pb
img0pb=os.path.basename(vis).replace('noise.aca.cycle6.','')+'.clean_'+'{1:s}'+'.pb'
img1pb=os.path.basename(vis).replace('noise.aca.cycle6.','')+'.clean_'+'{2:s}'+'.pb'
pbpdiff = img0+'.vs.'+img1.replace(shortvis+'.clean_','')+'_pbpdiff'
os.system('rm -rf '+pbpdiff)
immath(imagename=[img0pb,img1pb],
       expr='100*(IM0-IM1)/IM0',
       outfile=pbpdiff)

aU.imviewFields(imgs=[pbpdiff],
                contourImage=img0.replace('image.pbcor','pb'), 
                levels=[0.5], unit=1.0, contourThickness=2)

png1=img0+'.imview.png'
png2=img1+'.imview.png'
png3=diff+'.imview.png'
png4=pdiff+'.imview.png'
aU.montagePngs(png1=png1,png2=png2,png3=png3,png4=png4,tile='2x2',geometry='auto+10',
               outname=shortvis+'_clean_montage.png')

images=[img0,img1,diff,pdiff,img0pb,pbpdiff]

for i in images:    
    os.system('rm -rf '+i+'.fits')
    exportfits(i,i+'.fits')


print('')
print('The cleaned, pbcored annular rms noise for {1:s} vs {2:s} is: %.5f and %.5f'%(rms0,rms1))

print('The min, max, and median percent differences inside the 0.5 PB ({1:s}) level are: %.2f, %.2f, %.2f'%(statspdiff['min'][0],statspdiff['max'][0],statspdiff['median'][0]))

print('The min, max, and median percent differences between the 0.2 and  0.5 PB ({1:s}) level are: %.2f, %.2f, %.2f'%(outstatspdiff['min'][0],outstatspdiff['max'][0],outstatspdiff['median'][0]))

print('')

print('The min, max, and median differences inside the 0.5 PB ({1:s}) level are: %.3f, %.3f, %.3f'%(statsdiff['min'][0],statsdiff['max'][0],statsdiff['median'][0]))

print('The min, max, and median differences between the 0.2 and  0.5 PB ({1:s}) level are: %.3f, %.3f, %.3f'%(outstatsdiff['min'][0],outstatsdiff['max'][0],outstatsdiff['median'][0]))

with open('output.csv','w') as csvfile:    
    mywriter = csv.writer(csvfile)
    mywriter.writerow(['# vis',
                       'casa1_name', 'casa1_path', 
                       'casa2_name', 'casa2_path',
                       'x_size', 'y_size',
                       'array',
                       'rms_{1:s}', 'rms_{2:s}',
                       'pdiff_min_pb05', 'pdiff_max_pb05', 'pdiff_med_pb05',
                       'pdiff_min_pb0502', 'pdiff_max_pb0502', 'pdiff_med_pb0502',
                       'diff_min_pb05', 'diff_max_pb05', 'diff_med_pb05',
                       'diff_min_pb0502', 'diff_max_pb0502', 'diff_med_pb0502'])
    mywriter.writerow(['{0:s}', 
                       '{1:s}', '{3:s}',
                       '{2:s}', '{4:s}',
                       '{5:d}','{6:d}',
                       '{7:d}',
                       rms0,rms1,
                       statspdiff['min'][0],statspdiff['max'][0],statspdiff['median'][0],
                       outstatspdiff['min'][0],outstatspdiff['max'][0],outstatspdiff['median'][0],
                       statsdiff['min'][0],statsdiff['max'][0],statsdiff['median'][0],
                       outstatsdiff['min'][0],outstatsdiff['max'][0],outstatsdiff['median'][0]])


    '''.format(vis, casa1['name'], casa2['name'],casa1['path'],casa2['path'],simsize[0],simsize[1],array)
    
    fout = open('analysis.py','w')
    fout.write('#CASA1: '+casa1['path']+'\n')
    fout.write('#CASA2: '+casa1['path']+'\n')
    fout.write(outStr)
    fout.close()

def collateResults(dirlist,outcsv):
    '''
   
    Purpose: collate results from individual tests into single csv file
   
    '''

    import re

    firsthdr = True

    fileout = open(outcsv,'w')

    for mydir in dirlist:
        filename = os.path.join(mydir,'output.csv')
        if os.path.exists(filename):
            filein = open(filename,'r')
            for line in filein:
                if re.match('#',line):
                    if firsthdr:
                        fileout.write(line)
                        firsthdr = False
                    else:
                        continue
                else:
                    fileout.write(line)                                
            filein.close()
    fileout.close()



