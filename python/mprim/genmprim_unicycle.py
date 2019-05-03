# % /*
# %  * Copyright (c) 2008, Maxim Likhachev
# %  * All rights reserved.
# %  *
# %  * Redistribution and use in source and binary forms, with or without
# %  * modification, are permitted provided that the following conditions are met:
# %  *
# %  *     * Redistributions of source code must retain the above copyright
# %  *       notice, this list of conditions and the following disclaimer.
# %  *     * Redistributions in binary form must reproduce the above copyright
# %  *       notice, this list of conditions and the following disclaimer in the
# %  *       documentation and/or other materials provided with the distribution.
# %  *     * Neither the name of the Carnegie Mellon University nor the names of its
# %  *       contributors may be used to endorse or promote products derived from
# %  *       this software without specific prior written permission.
# %  *
# %  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# %  * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# %  * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# %  * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# %  * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# %  * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# %  * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# %  * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# %  * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# %  * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# %  * POSSIBILITY OF SUCH DAMAGE.
# %  */
#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

def genmprim_unicycle(outfilename):
# %
# %generates motion primitives and saves them into file
# %
# %written by Maxim Likhachev
# %---------------------------------------------------
# %

# %defines
    plt.figure()
    UNICYCLE_MPRIM_16DEGS = 1

    if UNICYCLE_MPRIM_16DEGS == 1:
        resolution = 0.025
        numberofangles = 16 #; %preferably a power of 2, definitely multiple of 8
        numberofprimsperangle = 5
        # %multipliers (multiplier is used as costmult*cost)
        forwardcostmult = 1
        backwardcostmult = 5
        forwardandturncostmult = 2
        sidestepcostmult = 10
        turninplacecostmult = 5
        # %note, what is shown x,y,theta changes (not absolute numbers)
        # %0 degreees
        basemprimendpts0_c = np.zeros((numberofprimsperangle, 4))#; %x,y,theta,costmult
        # %x aligned with the heading of the robot, angles are positive
        #%counterclockwise
        # %0 theta change
        basemprimendpts0_c[0,:] = [1, 0, 0, forwardcostmult]
        basemprimendpts0_c[1,:] = [8, 0, 0, forwardcostmult]
        basemprimendpts0_c[2,:] = [-1, 0, 0, backwardcostmult]
        # %1/16 theta change
        basemprimendpts0_c[3,:] = [8, 1, 1, forwardandturncostmult]
        basemprimendpts0_c[4,:] = [8, -1, -1, forwardandturncostmult]
        # %turn in place
        # %basemprimendpts0_c(6,:) = [0 0 1 turninplacecostmult];
        # %basemprimendpts0_c(7,:) = [0 0 -1 turninplacecostmult];
        # %45 degrees
        basemprimendpts45_c = np.zeros((numberofprimsperangle, 4))#; %x,y,theta,costmult (multiplier is used as costmult*cost)
        # %x aligned with the heading of the robot, angles are positiv
        # %counterclockwise
        # %0 theta change
        basemprimendpts45_c[0,:] = [1, 1, 0, forwardcostmult]
        basemprimendpts45_c[1,:] = [6, 6, 0, forwardcostmult]
        basemprimendpts45_c[2,:] = [-1, -1, 0, backwardcostmult]
        # %1/16 theta change
        basemprimendpts45_c[3,:] = [5, 7, 1, forwardandturncostmult]
        basemprimendpts45_c[4,:] = [7, 5, -1, forwardandturncostmult];
        # %turn in place
        # %basemprimendpts45_c(6,:) = [0 0 1 turninplacecostmult];
        #%basemprimendpts45_c(7,:) = [0 0 -1 turninplacecostmult];

        # %22.5 degrees
        basemprimendpts22p5_c = np.zeros((numberofprimsperangle, 4))#; %x,y,theta,costmult (multiplier is used as costmult*cost)
        # %x aligned with the heading of the robot, angles are positive
        # %counterclockwise
        # %0 theta change
        basemprimendpts22p5_c[0,:] = [2, 1, 0, forwardcostmult]
        basemprimendpts22p5_c[1,:] = [6, 3, 0, forwardcostmult]
        basemprimendpts22p5_c[2,:] = [-2, -1, 0, backwardcostmult]
        #%1/16 theta change
        basemprimendpts22p5_c[3,:] = [5, 4, 1, forwardandturncostmult]
        basemprimendpts22p5_c[4,:] = [7, 2, -1, forwardandturncostmult]
        # %turn in place
        # %basemprimendpts22p5_c(6,:) = [0 0 1 turninplacecostmult];
        # %basemprimendpts22p5_c(7,:) = [0 0 -1 turninplacecostmult];
    else:
        print('ERROR: undefined mprims type')
        return

    file = open(outfilename, 'w');
    #%write the header
    file.write('resolution_m: %f\n' % resolution)
    file.write('numberofangles: %d\n' %numberofangles)
    file.write('totalnumberofprimitives: %d\n' % numberofprimsperangle*numberofangles);

    # %iterate over angles
    for angleind in range(numberofangles):
        #figure(1);
        #hold off;
        #text(0, 0, int2str(angleind));
        #%iterate over primitives
        for primind in range(numberofprimsperangle):
            file.write('primID: %d \n' % primind)
            file.write('startangle_c: %d\n' % angleind);

        #%current angle
        currentangle = (angleind)*2*np.pi/numberofangles;
        currentangle_36000int = np.round((angleind)*36000/numberofangles);

        # %compute which template to use
        if (np.remainder(currentangle_36000int, 9000) == 0):
            basemprimendpts_c = basemprimendpts0_c[primind,:]
            angle = currentangle;
        elif (np.remainder(currentangle_36000int, 4500) == 0):
            basemprimendpts_c = basemprimendpts45_c[primind,:]
            angle = currentangle - 45*np.pi/180;
        elif (np.remainder(currentangle_36000int-7875, 9000) == 0):
            basemprimendpts_c = basemprimendpts33p75_c[primind,:]
            basemprimendpts_c[0] = basemprimendpts33p75_c[primind, 2]# %reverse x and y
            basemprimendpts_c[1] = basemprimendpts33p75_c[primind, 1]
            basemprimendpts_c[2] = -basemprimendpts33p75_c[primind, 3] # %reverse the angle as well
            angle = currentangle - 78.75*np.pi/180
        elif (np.remainder(currentangle_36000int-6750, 9000) == 0):
            basemprimendpts_c = basemprimendpts22p5_c[primind,:]
            basemprimendpts_c[0] = basemprimendpts22p5_c[primind, 1] #%reverse x and y
            basemprimendpts_c[1] = basemprimendpts22p5_c[primind, 0]
            basemprimendpts_c[2] = -basemprimendpts22p5_c[primind, 2] #%reverse the angle as well
            # %fprintf(1, '%d %d %d onto %d %d %d\n', basemprimendpts22p5_c(1), basemprimendpts22p5_c(2), basemprimendpts22p5_c(3), ...
            # %    basemprimendpts_c(1), basemprimendpts_c(2), basemprimendpts_c(3));
            angle = currentangle - 67.5*np.pi/180;
        elif (np.remainder(currentangle_36000int-5625, 9000) == 0):
            basemprimendpts_c = basemprimendpts11p25_c[primind,:]#;
            basemprimendpts_c[0] = basemprimendpts11p25_c[primind, 1]# %reverse x and y
            basemprimendpts_c[1] = basemprimendpts11p25_c[primind, 0]#;
            basemprimendpts_c[2] = -basemprimendpts11p25_c[primind, 2]#; %reverse the angle as well
            angle = currentangle - 56.25*np.pi/180;
        elif (np.remainder(currentangle_36000int-3375, 9000) == 0):
            basemprimendpts_c = basemprimendpts33p75_c[primind,:]
            angle = currentangle - 33.75*np.pi/180;
        elif (np.remainder(currentangle_36000int-2250, 9000) == 0):
            basemprimendpts_c = basemprimendpts22p5_c[primind,:]
            angle = currentangle - 22.5*np.pi/180;
        elif (np.remainder(currentangle_36000int-1125, 9000) == 0):
            basemprimendpts_c = basemprimendpts11p25_c[primind,:]
            angle = currentangle - 11.25*np.pi/180;
        else:
            print 'ERROR: invalid angular resolution. angle = %d' % currentangle_36000int
            return;

        #%now figure out what action will be
        baseendpose_c = basemprimendpts_c[0:3]
        additionalactioncostmult = basemprimendpts_c[3]
        endx_c = np.round(baseendpose_c[0]*np.cos(angle) - baseendpose_c[1]*np.sin(angle));
        endy_c = np.round(baseendpose_c[0]*np.sin(angle) + baseendpose_c[1]*np.cos(angle));
        endtheta_c = np.remainder(angleind + baseendpose_c[2], numberofangles);
        endpose_c = [endx_c, endy_c, endtheta_c]

        #fprintf(1, 'rotation angle=%f\n', angle*180/pi);

        #if baseendpose_c[1] == 0 and baseendpose_c[2] == 0
        #    %fprintf(1, 'endpose=%d %d %d\n', endpose_c(1), endpose_c(2), endpose_c(3));
        #end;

        #%generate intermediate poses (remember they are w.r.t 0,0 (and not
        #%centers of the cells)
        numofsamples = 10;
        intermcells_m = np.zeros((numofsamples,3))
        if UNICYCLE_MPRIM_16DEGS == 1:
            startpt = [0, 0 ,currentangle];
            endpt = [endpose_c[0]*resolution, \
                     endpose_c[1]*resolution, \
                     np.remainder(angleind + baseendpose_c[2], numberofangles)*2*np.pi/numberofangles]
            if ((endx_c == 0 and endy_c == 0) or baseendpose_c[2] == 0): # %turn in place or move forward
                for iind in range(0,numofsamples):
                    intermcells_m[iind,:] = [startpt[0] + (endpt[0] - startpt[0]*(iind)/(numofsamples-1)),
                                            startpt[1] + (endpt[1] - startpt[1]*(iind)/(numofsamples-1)),
                                            0]
                    rotation_angle = (baseendpose_c[2] ) * (2*np.pi/numberofangles);
                    intermcells_m[iind,2] = np.remainder(startpt[2] + (rotation_angle)*(iind)/(numofsamples-1), 2*np.pi);
            else: #%unicycle-based move forward or backward
                R = np.array([[np.cos(startpt[2]), np.sin(endpt[2]) - np.sin(startpt[2])],
                              [np.sin(startpt[2]), -np.cos(endpt[2]) - np.cos(startpt[2])]])
                S = np.dot(np.linalg.pinv(R),np.array([[endpt[0] - startpt[0]], [endpt[1] - startpt[1]]]))
                l = S[0];
                tvoverrv = S[1]
                rv = (baseendpose_c[2]*2*np.pi/numberofangles + l/tvoverrv);
                tv = tvoverrv*rv;

                if l < 0:
                    print 'WARNING: l = %d < 0 -> bad action start/end points' % l
                    l = 0;
                #%compute rv
                #%rv = baseendpose_c(3)*2*pi/numberofangles;
                #%compute tv
                #%tvx = (endpt(1) - startpt(1))*rv/(sin(endpt(3)) - sin(startpt(3)))
                #%tvy = -(endpt(2) - startpt(2))*rv/(cos(endpt(3)) - cos(startpt(3)))
                #%tv = (tvx + tvy)/2.0;
                #%generate samples
                for iind in range(0,numofsamples):
                    dt = (iind)/(numofsamples);

                    #%dtheta = rv*dt + startpt[2];
                    #%intermcells_m(iind,:) = [startpt(1) + tv/rv*(sin(dtheta) - sin(startpt(3))) ...
                    #%                        startpt(2) - tv/rv*(cos(dtheta) - cos(startpt(3))) ...
                    #%                        dtheta];

                    if(dt*tv < l):
                        intermcells_m[iind,:] = [startpt[0] + dt*tv*np.cos(startpt[2]),
                                                 startpt[1] + dt*tv*np.sin(startpt[2]),
                                                 startpt[2]];
                    else:
                        dtheta = rv*(dt - l/tv) + startpt[2];
                        intermcells_m[iind,:] = [startpt[0] + l*np.cos(startpt[2]) + tvoverrv*(np.sin(dtheta) - np.sin(startpt[2])),
                                                 startpt[1] + l*np.sin(startpt[2]) - tvoverrv*(np.cos(dtheta) - np.cos(startpt[2])),
                                                 dtheta];

                #%correct
                errorxy = [endpt[0] - intermcells_m[numofsamples-1,0] ,
                           endpt[1] - intermcells_m[numofsamples-1,1]]
                print  1.0/numofsamples
                #fprintf(1, 'l=%f errx=%f erry=%f\n', l, errorxy(1), errorxy(2));
                interpfactor = np.arange(0,1, 1.0/numofsamples)#[0:1/(numofsamples-1):1];
                intermcells_m[:,0] = intermcells_m[:,0] + errorxy[0]*interpfactor
                intermcells_m[:,1] = intermcells_m[:,1] + errorxy[1]*interpfactor


        #%write out

        file.write('endpose_c: %d %d %d\n' % (endpose_c[0], endpose_c[1], endpose_c[2]))
        #%write out
        file.write('additionalactioncostmult: %d\n'% additionalactioncostmult)
        file.write('intermediateposes: %d \n'% intermcells_m.shape[0])
        for interind in range (0, intermcells_m.shape[0]):
            print "i ", interind
            file.write('%.4f %.4f %.4f\n' % (intermcells_m[interind,0], intermcells_m[interind,1], intermcells_m[interind,2]))

        plt.plot(intermcells_m[:,0], intermcells_m[:,1])
        plt.show(False)
        plt.waitforbuttonpress()
    file.close()


if __name__ == '__main__':
    genmprim_unicycle("my_unicycle.mprim")
