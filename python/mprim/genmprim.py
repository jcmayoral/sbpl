#! /usr/bin/python
import numpy as np
import matplotlib.pyplot as plt

def genmprim(outfilename):

    # %
    # %generates motion primitives and saves them into file
    # %
    # %written by Maxim Likhachev
    # %---------------------------------------------------
    # %
    #
    # %defines
    plt.figure()

    LINESEGMENT_MPRIMS = 1 #%set the desired type of motion primitives
    UNICYCLE_MPRIMS = 0


    if LINESEGMENT_MPRIMS == 1:
        resolution = 0.01
        numberofangles = 16#32#; %preferably a power of 2, definitely multiple of 8
        numberofprimsperangle = 16#;

        #%multipliers (multiplier is used as costmult*cost)
        forwardcostmult = 1
        backwardcostmult = 5
        forwardandturncostmult = 1
        sidestepcostmult = 50
        turninplacecostmult = 50

        # %note, what is shown x,y,theta changes (not absolute numbers)

        # %0 degreees
        basemprimendpts0_c = np.zeros((numberofprimsperangle, 4))#; %x,y,theta,costmult
        # %x aligned with the heading of the robot, angles are positive
        # %counterclockwise
        # %0 theta change
        basemprimendpts0_c[0,:] = [2,0,0,forwardcostmult]
        basemprimendpts0_c[1,:] = [8,0,0,forwardcostmult]
        basemprimendpts0_c[2,:] = [16,0,0,forwardcostmult]
        basemprimendpts0_c[3,:] = [12,2,0,sidestepcostmult]
        basemprimendpts0_c[4,:] = [12,-2,0,sidestepcostmult]
        basemprimendpts0_c[5,:] = [4,3,0,sidestepcostmult]
        basemprimendpts0_c[6,:] = [4,-3,0,sidestepcostmult]
        basemprimendpts0_c[7,:] = [-10,0,0,backwardcostmult]
        # %1/32 theta change
        basemprimendpts0_c[8,:] = [12,2,1,forwardandturncostmult]
        basemprimendpts0_c[9,:] = [12,-2,-1,forwardandturncostmult]
        # %2/32 theta change
        basemprimendpts0_c[10,:] = [8,3,2,forwardandturncostmult]
        basemprimendpts0_c[11,:] = [8,-3,-2,forwardandturncostmult]
        # %turn in place
        basemprimendpts0_c[12,:] = [0,0,1,turninplacecostmult]
        basemprimendpts0_c[13,:] = [0,0,-1,turninplacecostmult]
        basemprimendpts0_c[14,:] = [0,0,3,turninplacecostmult]
        basemprimendpts0_c[15,:] = [0,0,-3,turninplacecostmult]

        # %45 degrees
        basemprimendpts45_c = np.zeros((numberofprimsperangle, 4))#; %x,y,theta,costmult (multiplier is used as costmult*cost)
        # %x aligned with the heading of the robot, angles are positive
        # %counterclockwise
        # %0 theta change
        basemprimendpts45_c[0,:] = [2,1,0,forwardcostmult]
        basemprimendpts45_c[1,:] = [6,3,0,forwardcostmult]
        basemprimendpts45_c[2,:] = [12,6,0,forwardcostmult]
        basemprimendpts45_c[3,:] = [4,6,0,sidestepcostmult]
        basemprimendpts45_c[4,:] = [12,2,0,sidestepcostmult]
        basemprimendpts45_c[5,:] = [0,4,0,sidestepcostmult]
        basemprimendpts45_c[6,:] = [8,0,0,sidestepcostmult]
        basemprimendpts45_c[7,:] = [-8,-4,0,backwardcostmult]
        # %1/32 theta change
        basemprimendpts45_c[8,:] = [4,6,1,forwardandturncostmult]
        basemprimendpts45_c[9,:] = [12,2,-1,forwardandturncostmult]
        # %2/32 theta change
        basemprimendpts45_c[10,:] = [2,5,2,forwardandturncostmult]
        basemprimendpts45_c[11,:] = [10,1,-2,forwardandturncostmult]
        # %turn in place
        basemprimendpts45_c[12,:] = [0,0,1,turninplacecostmult]
        basemprimendpts45_c[13,:] = [0,0,-1,turninplacecostmult]
        basemprimendpts45_c[14,:] = [0,0,3,turninplacecostmult]
        basemprimendpts45_c[15,:] = [0,0,-3,turninplacecostmult]

        # %22.5 degrees
        basemprimendpts22p5_c = np.zeros((numberofprimsperangle, 4))#; %x,y,theta,costmult (multiplier is used as costmult*cost)
        # %x aligned with the heading of the robot, angles are positive
        # %counterclockwise
        # %0 theta change
        basemprimendpts22p5_c[0,:] = [4,1,0,forwardcostmult]
        basemprimendpts22p5_c[1,:] = [8,2,0,forwardcostmult]
        basemprimendpts22p5_c[2,:] = [23,3,0,forwardcostmult]
        basemprimendpts22p5_c[3,:] = [8,4,0,sidestepcostmult]
        basemprimendpts22p5_c[4,:] = [23,2,0,sidestepcostmult]
        basemprimendpts22p5_c[5,:] = [0,3,0,sidestepcostmult]
        basemprimendpts22p5_c[6,:] = [8,-1,0,sidestepcostmult]
        basemprimendpts22p5_c[7,:] = [-8,-2,0,backwardcostmult]
        # %1/32 theta change
        basemprimendpts22p5_c[8,:] = [8,4,1,forwardandturncostmult]
        basemprimendpts22p5_c[9,:] = [12,2,-1,forwardandturncostmult]
        # %2/32 theta change
        basemprimendpts22p5_c[10,:] = [4,4,2,forwardandturncostmult]
        basemprimendpts22p5_c[11,:] = [12,0,-2,forwardandturncostmult]
        # %turn in place
        basemprimendpts22p5_c[12,:] = [0,0,1,turninplacecostmult]
        basemprimendpts22p5_c[13,:] = [0,0,-1,turninplacecostmult]
        basemprimendpts22p5_c[14,:] = [0,0,3,turninplacecostmult]
        basemprimendpts22p5_c[15,:] = [0,0,-3,turninplacecostmult]

        # %11.25 degrees
        basemprimendpts11p25_c = np.zeros((numberofprimsperangle, 4))#; %x,y,theta,costmult (multiplier is used as costmult*cost)
        # %x aligned with the heading of the robot, angles are positive
        # %counterclockwise
        # %0 theta change
        basemprimendpts11p25_c[0,:] = [6,1,0,forwardcostmult]
        basemprimendpts11p25_c[1,:] = [12,2,0,forwardcostmult]
        basemprimendpts11p25_c[2,:] = [18,3,0,forwardcostmult]
        basemprimendpts11p25_c[3,:] = [8,3,0,sidestepcostmult]
        basemprimendpts11p25_c[4,:] = [12,0,0,sidestepcostmult]
        basemprimendpts11p25_c[5,:] = [2,3,0,sidestepcostmult]
        basemprimendpts11p25_c[6,:] = [6,-2,0,sidestepcostmult]
        basemprimendpts11p25_c[7,:] = [-12,-2,0,backwardcostmult]
        # %1/32 theta change
        basemprimendpts11p25_c[8,:] = [8,3,1,forwardandturncostmult]
        basemprimendpts11p25_c[9,:] = [12,0,-1,forwardandturncostmult]
        # %2/32 theta change
        basemprimendpts11p25_c[10,:] = [4,4,2,forwardandturncostmult]
        basemprimendpts11p25_c[11,:] = [10,-1,-2,forwardandturncostmult]
        # %turn in place
        basemprimendpts11p25_c[12,:] = [0,0,1,turninplacecostmult]
        basemprimendpts11p25_c[13,:] = [0,0,-1,turninplacecostmult]
        basemprimendpts11p25_c[14,:] = [0,0,3,turninplacecostmult]
        basemprimendpts11p25_c[15,:] = [0,0,-3,turninplacecostmult]

        # %33.75 degrees
        basemprimendpts33p75_c = np.zeros((numberofprimsperangle, 4)) # %x,y,theta,costmult
        # %x aligned with the heading of the robot, angles are positive
        # %counterclockwise
        # %0 theta change
        basemprimendpts33p75_c[0,:] = [6,2,0,forwardcostmult]
        basemprimendpts33p75_c[1,:] = [12,4,0,forwardcostmult]
        basemprimendpts33p75_c[2,:] = [18,6,0,forwardcostmult]
        basemprimendpts33p75_c[3,:] = [8,5,0,sidestepcostmult]
        basemprimendpts33p75_c[4,:] = [12,2,0,sidestepcostmult]
        basemprimendpts33p75_c[5,:] = [0,4,0,sidestepcostmult]
        basemprimendpts33p75_c[6,:] = [6,-2,0,sidestepcostmult]
        basemprimendpts33p75_c[7,:] = [-12,-4,0,backwardcostmult]
        # %1/32 theta change
        basemprimendpts33p75_c[8,:] = [8,5,1,forwardandturncostmult]
        basemprimendpts33p75_c[9,:] = [12,2,-1,forwardandturncostmult]
        # %2/32 theta change
        basemprimendpts33p75_c[10,:] = [2,5,2,forwardandturncostmult]
        basemprimendpts33p75_c[11,:] = [6,-2,-2,forwardandturncostmult]
        # %turn in place
        basemprimendpts33p75_c[12,:] = [0,0,1,turninplacecostmult]
        basemprimendpts33p75_c[13,:] = [0,0,-1,turninplacecostmult]
        basemprimendpts33p75_c[14,:] = [0,0,3,turninplacecostmult]
        basemprimendpts33p75_c[15,:] = [0,0,-3,turninplacecostmult]


    elif UNICYCLE_MPRIMS == 1:
        print 'ERROR: unsupported mprims type'
        return
    else:
        print 'ERROR: undefined mprims type'
        return

    fout = open(outfilename, 'w');
    #%write the header
    fout.write('resolution_m: %f\n' % resolution)
    fout.write('numberofangles: %d\n' % numberofangles)
    fout.write('totalnumberofprimitives: %d\n' % int(numberofprimsperangle*numberofangles))

    #%iterate over angles
    print numberofangles
    for angleind in range(numberofangles):
        print angleind

        #figure(1);
        #hold off;
        #text(0, 0, int2str(angleind));

        #%iterate over primitives
        for primind in range(numberofprimsperangle):
            fout.write('primID: %d \n'% int(primind))
            fout.write('startangle_c: %d\n'% int(angleind))

            #%current angle
            currentangle = (angleind)*2*np.pi/numberofangles;
            currentangle_36000int = np.round((angleind)*36000/numberofangles);

            #%compute which template to use
            if (np.remainder(currentangle_36000int, 9000) == 0):
                basemprimendpts_c = basemprimendpts0_c[primind,:]
                angle = currentangle
            elif (np.remainder(currentangle_36000int, 4500) == 0):
                basemprimendpts_c = basemprimendpts45_c[primind,:]
                angle = currentangle - 45*np.pi/180;
            elif (np.remainder(currentangle_36000int-7875, 9000) == 0):
                basemprimendpts_c = basemprimendpts33p75_c[primind,:]
                basemprimendpts_c[0] = basemprimendpts33p75_c[primind, 1] #; %reverse x and y
                basemprimendpts_c[1] = basemprimendpts33p75_c[primind, 0]#;
                basemprimendpts_c[2] = -basemprimendpts33p75_c[primind, 3]#"; %reverse the angle as well
                angle = currentangle - 78.75*np.pi/180
            elif (np.remainder(currentangle_36000int-6750, 9000) == 0):
                basemprimendpts_c = basemprimendpts22p5_c[primind,:]
                basemprimendpts_c[0] = basemprimendpts22p5_c[primind, 1]#; %reverse x and y
                basemprimendpts_c[1] = basemprimendpts22p5_c[primind, 0]#;
                basemprimendpts_c[2] = -basemprimendpts22p5_c[primind, 2] #%reverse the angle as well
                #%fprintf(1, '%d %d %d onto %d %d %d\n', basemprimendpts22p5_c(1), basemprimendpts22p5_c(2), basemprimendpts22p5_c(3), ...
                #%    basemprimendpts_c(1), basemprimendpts_c(2), basemprimendpts_c(3));
                angle = currentangle - 67.5*np.pi/180;
            elif (np.remainder(currentangle_36000int-5625, 9000) == 0):
                basemprimendpts_c = basemprimendpts11p25_c[primind,:]
                basemprimendpts_c[0] = basemprimendpts11p25_c[primind, 1] #%reverse x and y
                basemprimendpts_c[1] = basemprimendpts11p25_c[primind, 0]
                basemprimendpts_c[2] = -basemprimendpts11p25_c[primind, 2]#; %reverse the angle as well
                angle = currentangle - 56.25*np.pi/180
            elif (np.remainder(currentangle_36000int-3375, 9000) == 0):
                basemprimendpts_c = basemprimendpts33p75_c[primind,:]
                angle = currentangle - 33.75*np.pi/180
            elif (np.remainder(currentangle_36000int-2250, 9000) == 0):
                basemprimendpts_c = basemprimendpts22p5_c[primind,:]
                angle = currentangle - 22.5*np.pi/180
            elif (np.remainder(currentangle_36000int-1125, 9000) == 0):
                basemprimendpts_c = basemprimendpts11p25_c[primind,:]
                angle = currentangle - 11.25*np.pi/180;
            else:
                print 'ERROR: invalid angular resolution. angle = %d' % currentangle_36000int
                return;

            # %now figure out what action will be
            baseendpose_c = basemprimendpts_c[0:3]
            additionalactioncostmult = basemprimendpts_c[3]
            endx_c = np.round(baseendpose_c[0]*np.cos(angle) - baseendpose_c[1]*np.sin(angle))
            endy_c = np.round(baseendpose_c[1]*np.sin(angle) + baseendpose_c[1]*np.cos(angle))
            endtheta_c = np.remainder(angleind + baseendpose_c[2], numberofangles);
            endpose_c = [endx_c,endy_c,endtheta_c];


            # if baseendpose_c(2) == 0 & baseendpose_c(3) == 0
                # %fprintf(1, 'endpose=%d %d %d\n', endpose_c(1), endpose_c(2), endpose_c(3));
            # end;

            # %generate intermediate poses (remember they are w.r.t 0,0 (and not
            # %centers of the cells)
            numofsamples = 10;
            intermcells_m = np.zeros((numofsamples,3))
            if LINESEGMENT_MPRIMS == 1:
                startpt = [0, 0, currentangle];
                endpt = [endpose_c[0]*resolution,  \
                        endpose_c[1]*resolution,   \
                        np.remainder(angleind + baseendpose_c[2], numberofangles)*2*np.pi/numberofangles]
                intermcells_m = np.zeros((numofsamples,3));
                for iind in range(0,numofsamples):
                    intermcells_m[iind,:] = [startpt[0] + (endpt[0] - startpt[0])*(iind)/(numofsamples-1),
                                            startpt[1] + (endpt[1] - startpt[1])*(iind)/(numofsamples-1),
                                            0]
                    rotation_angle = (baseendpose_c[2] ) * (2*np.pi/numberofangles)
                    intermcells_m[iind,2] = np.remainder(startpt[2] + (rotation_angle)*(iind)/(numofsamples-1), 2*np.pi)

            fout.write('endpose_c: %d %d %d\n' % (endpose_c[0], endpose_c[1], endpose_c[2]))
            #%write out
            fout.write('additionalactioncostmult: %d\n'% additionalactioncostmult)
            print intermcells_m.shape[0]
            fout.write('intermediateposes: %d \n'% intermcells_m.shape[0])
            for interind in range (0, intermcells_m.shape[0]):
                print "i ", interind
                fout.write('%.4f %.4f %.4f\n' % (intermcells_m[interind,0], intermcells_m[interind,1], intermcells_m[interind,2]))

            plt.plot(intermcells_m[:,0], intermcells_m[:,1]);
            plt.show(False)
            plt.waitforbuttonpress()

            #text(intermcells_m(numofsamples,1), intermcells_m(numofsamples,2), int2str(endpose_c(3)));
            #hold on;

    fout.close()
    #fclose('all');


if __name__ == '__main__':
    genmprim("my_prim_test.mprim")
