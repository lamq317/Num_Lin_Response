function SolveSchro(Pot,pts,deltaq,number)

#Ham=zeros(ComplexF64,pts,pts);
Ham=zeros(Float64,pts,pts)

#Ham[1,1]=-2.0
#Ham[1,2]=1.0
Ham[1,1]=0.0
Ham[1,2]=0.0
    #Ham[0,pts-1]=1.0

#Ham[pts,pts]=-2.0
#Ham[pts,pts-1]=1.0
Ham[pts,pts]=0.0
Ham[pts,pts-1]=0.0

for i in 2:pts-1
    Ham[i,i]=-2.0
    Ham[i,i+1]=1.0
    Ham[i,i-1]=1.0

end

for i in 1:pts
    Ham[i,:]=-(1/(Pot[i]*deltaq^2))*Ham[i,:]
end

eigenValues=eigvals(Ham);
eigenVectors=eigvecs(Ham);


for i in 3:(number+2)
    RNorm=deltaq*dot(Pot.*eigenVectors[:,i],eigenVectors[:,i])
    RNorm=sqrt(RNorm)
    eigenVectors[:,i]=eigenVectors[:,i]/RNorm
end

    eigenValues[3:(number+2)],eigenVectors[:,3:(number+2)]

end


function NumGreensInc(grid,epsilon,psi, phi, k)

    kabs=abs(k)
    deltaq=abs(grid[2]-grid[1])

    gridpts=size(grid,1)
    #cdef int half=int(gridpts/2)

    integ1=convert(ComplexF64,0.0)
    #integ2=convert(ComplexF64,0.0)

    #cdef complex psi=np.conjugate(psi)

    for i in 1:gridpts

        #for j in 1:i-1
         for j in 1:gridpts
            integ1+=epsilon[i]*conj(psi[i])*exp(im*kabs*abs(grid[i]-grid[j]))*phi[j] #changed, sign of exponent

            #dum1=epsilon[i]*psi[i]*np.cos(k*(grid[i]-grid[j]))*phi[j]
            #integ1+=dum1

            #integ2+=dum1
            #dum1=epsilon[i]*conj(psi[i])*cos(k*(grid[i]-grid[j]))*phi[j]
            #dum2=epsilon[j]*conj(psi[j])*cos(k*(grid[j]-grid[i]))*phi[i]
            #integ1+=dum1-dum2
            #integ2+=dum1+dum2

        end
    end

    #adding the missing terms
    #for i in 1:gridpts
    #    integ2+=epsilon[i]*conj(psi[i])*phi[i]
    #end

    integ1=-integ1*deltaq^2 *im/kabs

    integ1

    #integ2=-2 *im*integ2*deltaq^2

    #integ1+integ2


end

function MaxwellInner(grid,epsilon,psi,phi)
    integrand=convert(ComplexF64,0.0)
    gridpts=size(grid,1)
    deltaq=abs(grid[1]-grid[2])

    for i in 1:gridpts
        integrand+=epsilon[i]*conj(psi[i])*phi[i]
    end

   integrand*deltaq
end

#auxiliary function to compute the inverse of free-partcile greens function
#SystMod is an array of size gridpts x dim. gridpts is size of grid.
#WORKING, but with some small deviations from python code...
function InvGreens(grid,epsilon,SystMod,k)

    gridpts,dim=size(SystMod)
    MatA=zeros(ComplexF64,(dim,dim))

    for i in 1:dim
        for j in 1:dim
            MatA[i,j]=NumGreensInc(grid,epsilon,SystMod[:,i],SystMod[:,j],k)
        end
    end


    inv(MatA)
end


#VERIFIED
function ActRetGreensPos(grid,phi,rfs,k)

    kabs=abs(k)
    deltaq=abs(grid[2]-grid[1])

    gridpts=size(grid,1)
    pts_rfs=size(rfs,1)



    integ1=convert(ComplexF64,0.0)
    #integ2=convert(ComplexF64,0.0)

    final=zeros(ComplexF64,pts_rfs)

    for i in 1:pts_rfs

        for j in 1:gridpts
            integ1+=exp(im*kabs*abs(rfs[i]-grid[j]))*phi[j] #changed sign of exponent

        end



       final[i]=-im*integ1*deltaq/kabs
       integ1=0.0




    end

    final


end

function TorthoM(grid,epsilon,SystMods,k)
    kabs=abs(k)
    gridpts,dim=size(SystMods)
    freeK=(1/sqrt(2*pi*kabs))*exp.(im*kabs*grid) #We are using Energy normalized free particle states!!
    freeminusK=(1/sqrt(2*pi*kabs))*exp.(-im*kabs*grid)

    Snk=zeros(ComplexF64,dim)
    SMat=zeros(ComplexF64,(2,2))

    #These ones are for the calculation of m=R modes
    Snkminus=zeros(ComplexF64,dim)

    for i in 1:dim
        Snk[i]=MaxwellInner(grid,epsilon,SystMods[:,i],freeK)


        Snkminus[i]=MaxwellInner(grid,epsilon,SystMods[:,i],freeminusK)
    end

    InvA=InvGreens(grid,epsilon,SystMods,k)

    #we use the ordering [[++,+-],[-+,--]]. dot dunction automatically takes the complex conjgate!!
    #SMat[1,1]=1+2*pi*im*dot(Snk,*(InvA,Snk))
    #SMat[1,2]=2*pi*im*dot(Snk,*(InvA,Snkminus))
    #SMat[2,1]=2*pi*im*dot(Snkminus,*(InvA,Snk))
    #SMat[2,2]=1+2*pi*im*dot(Snkminus,*(InvA,Snkminus))

    SMat[1,1]=-dot(Snk,*(InvA,Snk))
    SMat[1,2]=-dot(Snk,*(InvA,Snkminus))
    SMat[2,1]=-dot(Snkminus,*(InvA,Snk))
    SMat[2,2]=-dot(Snkminus,*(InvA,Snkminus))

    SMat



end


#This function returns a size(rfs) \times size(rfs) square matrix whose ijth entry corresponds to the
# 'direct scattering' resolvent evaluated at rfs[i],rfs[j]
#In this formulation, based on the position basis for the Maxwell potential, The rfs is an array
#that contains the positions that define  the potential
function DirectResolvPos(grid,epsilon,SystMods,rfs,InvA,k)
    kabs=abs(k)
    gridpts,dim=size(SystMods)
    rfs_pts=size(rfs,1)
    #SysG0phi=zeros(ComplexF64,(dim,1))
    result=zeros(ComplexF64,(rfs_pts,rfs_pts))
    KetContrib=zeros(ComplexF64,(rfs_pts,dim))
    BraContrib=zeros(ComplexF64,(dim,rfs_pts))

    #G0phi=ActRetGreens(grid,phi,k) #action of the retarded Green's function on phi

    #For each system mode, we compute its propagation by applying the free-particle retarded Greens function
    #store them in the Ket and Bra contributions
    for i in 1:dim

        KetContrib[:,i]=ActRetGreensPos(grid,SystMods[:,i],rfs,k)
        BraContrib[i,:]=ActRetGreensPos(grid,conj(SystMods[:,i]),rfs,k)
      # SysG0phi[i]=MaxwellInner(grid,epsilon,SystMods[:,i],G0phi)

    end

    for i in 1:rfs_pts

        for j in 1:rfs_pts
            result[i,j]=-im*exp(im*kabs*abs(rfs[i]-rfs[j]))/kabs #changed sign of exponent

        end


    end

    result+=-*(KetContrib,*(InvA,BraContrib))


   result

end


#Calculate the position representation of the inverse of Matrix M
#It is assumed that the array rfs contains all the spatial sampling of the potential V. Therefore, they must
#have the same size. IMPORTANT: the potential must be weighted by an extra epsilon factor, to account for
#the Maxwell inner product

function NonOverInvM(grid,epsilon,SystMods,rfs,InvA,V,k)
    gridpts,dim=size(SystMods)
    rfs_pts=size(rfs,1)
    M=zeros(ComplexF64,(rfs_pts,rfs_pts))

    G0rij=DirectResolvPos(grid,epsilon,SystMods,rfs,InvA,k)


    for i in 1:rfs_pts
        M[i,i]+=V[i] #as pointed out above, we assume the V array is weigthed by \epsilon

        for j in 1:rfs_pts
            M[i,j]+=-V[i]*G0rij[i,j]*V[j]
        end

    end

    inv(M) #


end




#Modified function to compute the value of KPlus at the space spanned by the potential.

#It returns the orthogonalized functions Kplus/Kminus evaluated at points rfs
#To avoid redundant calculations, this function returns 4 arrays, which read, from left to right:
#1. Kplus evaluated at rfs for |k|
#2. Kplus evaluated at rfs for -|k|
#3. Kminus evaluated at rfs for |k|
#4. Kminus evaluated at rfs for -|k|


function ComputeKplusPos(grid,rfs,epsilon,dim,SystMod,k)

    gridpts=size(grid,1)
    rfs_pts=size(rfs,1)
    column1Plus=zeros(ComplexF64,dim)
    #column2Plus=zeros(ComplexF64,dim)
    column1Minus=zeros(ComplexF64,dim)
    #column2Minus=zeros(ComplexF64,dim)


    resultPlusInc=zeros(ComplexF64,rfs_pts)
    resultMinusInc=zeros(ComplexF64,rfs_pts)
    resultPlusOut=zeros(ComplexF64,rfs_pts)
    resultMinusOut=zeros(ComplexF64,rfs_pts)

    kabs=abs(k)

    #generate real-space representation of free-particle with momentum k
    freepartPlus=sqrt(1/(2*pi*kabs))*exp.(im*kabs*grid)
    freepart_redPlus=sqrt(1/(2*pi*kabs))*exp.(im*kabs*rfs)

    freepartMinus=sqrt(1/(2*pi*kabs))*exp.(-im*kabs*grid)
    freepart_redMinus=sqrt(1/(2*pi*kabs))*exp.(-im*kabs*rfs)


    for i in 1:dim
        column1Plus[i]=MaxwellInner(grid,epsilon,SystMod[:,i],freepartPlus)
        column1Minus[i]=MaxwellInner(grid,epsilon,SystMod[:,i],freepartMinus)
        #row[i]=ActRetGreens(grid,SystMod[:,i],k)
    end

    InvAInc=InvGreens(grid,epsilon,SystMod,k)
    InvAOut=conj(InvAInc)

    columnPlusInc=*(InvAInc,column1Plus)
    columnMinusInc=*(InvAInc,column1Minus)
    columnPlusOut=*(InvAOut,column1Plus)
    columnMinusOut=*(InvAOut,column1Minus)

    #generate the final function
    for i in 1:dim
        #result+=column2[i]*ActRetGreens(grid,SystMod[:,i],k)
        resultPlusInc+=columnPlusInc[i]*ActRetGreensPos(grid,SystMods[:,i],rfs,k)
        resultMinusInc+=columnMinusInc[i]*ActRetGreensPos(grid,SystMods[:,i],rfs,k)
        resultPlusOut+=columnPlusOut[i]*conj(ActRetGreensPos(grid,SystMods[:,i],rfs,k))
        resultMinusOut+=columnMinusOut[i]*conj(ActRetGreensPos(grid,SystMods[:,i],rfs,k))
    end



    +(freepart_redPlus,-resultPlusInc),+(freepart_redMinus,-resultMinusInc),+(freepart_redPlus,-resultPlusOut),+(freepart_redMinus,-resultMinusOut)
end

#compute the background T matrix contribution to the background scattering matrix
function ComputeSBG(grid,rfs,epsilon,dim,SystMod,V,k)
    gridpts=size(grid,1)
    rfs_pts=size(rfs,1)
    Identity=zeros(ComplexF64,(2,2))
    Identity[1,1]=1
    Identity[2,2]=1

    ColPlusInc=zeros(ComplexF64,(rfs_pts,1))
    ColMinusInc=zeros(ComplexF64,(rfs_pts,1))
    RowPlusOut=zeros(ComplexF64,(1,rfs_pts))
    RowMinusOut=zeros(ComplexF64,(1,rfs_pts))

    TF=zeros(ComplexF64,(2,2))

    #calculation of the orthogonalized modes to the system modes at the points of the potential
    PlusInc,MinusInc,PlusOut,MinusOut=ComputeKplusPos(grid,rfs,epsilon,dim,SystMod,k)

    for i in 1:rfs_pts  #IMPORTANT: It is assumed that the potential already includes the permittivity and deltaq prefactors
        ColPlusInc[i,1]=V[i]*PlusInc[i]
        ColMinusInc[i,1]=V[i]*MinusInc[i]
        RowPlusOut[1,i]=V[i]*conj(PlusOut[i])
        RowMinusOut[1,i]=V[i]*conj(MinusOut[i])
    end

    InvA=InvGreens(grid,epsilon,SystMod,k)
    InvM=NonOverInvM(grid,epsilon,SystMod,rfs,InvA,V,k)

    #construction of the T matrix. We use the ordering [[++,+-],[-+,--]]
    TF[1,1]=*(RowPlusOut,*(InvM,ColPlusInc))[1,1]
    TF[1,2]=*(RowPlusOut,*(InvM,ColMinusInc))[1,1]
    TF[2,1]=*(RowMinusOut,*(InvM,ColPlusInc))[1,1]
    TF[2,2]=*(RowMinusOut,*(InvM,ColMinusInc))[1,1]

    #adding the 'orthogonalization' contribution
    TOrt=TorthoM(grid,epsilon,SystMod,k)
    #TOrt=zeros(ComplexF64,(2,2))
    Identity-2*pi*im*(TF+TOrt)

    #This the background scattering matrix
    #TBG=Identity-2*pi*im*(TOrt+TF)
    #TF


end


function calcSysBathNO(grid,rfs,V,epsilon,SystMods,ModFreqs,k)

    kabs=abs(k)
    gridpts,dim=size(SystMods)
    NsystDim=size(rfs,1)


    freeK=(1/sqrt(2*pi*kabs))*exp.(im*kabs*grid) #Using energy-normalized free-partcle states!
    freeminusK=(1/sqrt(2*pi*kabs))*exp.(-im*kabs*grid)
    #Building of all relevant vectors and matrices. Notice this is done such that we avoid redundant
    #calculations as much as possible:

    #vectors...
    Snk=zeros(ComplexF64,dim)
    Snkminus=zeros(ComplexF64,dim)

    VNsysKPlus=zeros(ComplexF64,NsystDim)
    VNsysKMinus=zeros(ComplexF64,NsystDim)



    #Matrices...
    #BMat=zeros(ComplexF64,(dim,dim))
    CMat=zeros(ComplexF64,(NsystDim,dim))
    CMatOut=zeros(ComplexF64,(NsystDim,dim)) #For the calculation of background scatetring matrix
    DMat=zeros(ComplexF64,(dim,NsystDim))

    Tbg=zeros(ComplexF64,(2,2))
    Identity=zeros(ComplexF64,(2,2))
    Identity[1,1]=1.0
    Identity[2,2]=1.0

    #To store final results
    Couplings=zeros(ComplexF64,(2,dim))


    for i in 1:dim
        Snk[i]=MaxwellInner(grid,epsilon,SystMods[:,i],freeK)

        Snkminus[i]=MaxwellInner(grid,epsilon,SystMods[:,i],freeminusK)

    end


    for i in 1:NsystDim

        VNsysKPlus[i]=V[i]*(1/sqrt(2*pi*kabs))*exp(im*kabs*rfs[i])
        VNsysKMinus[i]=V[i]*(1/sqrt(2*pi*kabs))*exp(-im*kabs*rfs[i])
    end

    for i in 1:dim

        dumG0Chi=ActRetGreensPos(grid,SystMods[:,i],rfs,k)
        dumG0Chiconj=conj(dumG0Chi)


        for l in 1:NsystDim

            CMat[l,i]=V[l]*dumG0Chi[l]
            CMatOut[l,i]=V[l]*dumG0Chiconj[l]
        end

    end


    for i in 1:dim

        dumG0Chi=ActRetGreensPos(grid,SystMods[:,i],rfs,k)

        for j in 1:NsystDim

            DMat[i,j]=dumG0Chi[j]*V[j]

        end

    end

    #calculation of InvA and InvM
    InvA=InvGreens(grid,epsilon,SystMods,k)
    InvAOut=conj(InvA)

    InvM=NonOverInvM(grid,epsilon,SystMods,rfs,InvA,V,k)


    #End of vector and matrix calculations *****

    #auxiliary vector to avoid redundant calculations
    VecPlusKBufInc=VNsysKPlus-*(CMat,*(InvA,Snk))            #useful for the calculation of the Tbg matrix
    VecPlusKminBufInc=VNsysKMinus-*(CMat,*(InvA,Snkminus))

    #For computation of background scattering matrix
    VecPlusKBufOut=VNsysKPlus-*(CMatOut,*(InvAOut,Snk))            #useful for the calculation of the Tbg matrix
    VecPlusKminBufOut=VNsysKMinus-*(CMatOut,*(InvAOut,Snkminus))

    VecPlusK=*(InvM,VecPlusKBufInc)
    VecPlusKmin=*(InvM,VecPlusKminBufInc)

    dumKin=*(InvA,*(DMat,VecPlusK))
    dumKinmin=*(InvA,*(DMat,VecPlusKmin))
    Kinetic=*(InvA,Snk)+dumKin
    KineticR=*(InvA,Snkminus)+dumKinmin



    #final step for calculation of coupling... DEPENDING ON THE REPRESENTATION OF BATH COORDINATES, WE WILL
                                                #HAVE A DIFFERENT SCALAR PREFACTOR

    Couplings[1,:]=Kinetic
    Couplings[2,:]=KineticR

    for i in 1:dim


        #Couplings[1,i]=Couplings[1,i]#/sqrt(ModFreqs[i])
        #Couplings[2,i]=Couplings[2,i]#/sqrt(ModFreqs[i])

        #VERY  IMPORTANT: note that throughout this code, I used the 'Scattering normalization' of the states
        #by diving the free particle states by \sqrt(2*pi*kabs). To get the coupling constants of the expressions developed
        #by Dominik and Viviescas, we need to divide by sqrt(\omega_{\lambda})

        Couplings[1,i]=Couplings[1,i]/sqrt(ModFreqs[i]) #This is \frac{\mathcal{W}_{\lambda,m}}{2 \sqrt{\omega \omega_{\lambda}}}
        Couplings[2,i]=Couplings[2,i]/sqrt(ModFreqs[i])





    end

    #calculation of Tbg matrix, we use the ordering [[++,+-],[-+,--]]. Remember that dot function
    #takes complex conjugate on left vector
    Tbg[1,1]=-dot(Snk,*(InvA,Snk))+dot(VecPlusKBufOut,*(InvM,VecPlusKBufInc))
    Tbg[1,2]=-dot(Snk,*(InvA,Snkminus))+dot(VecPlusKBufOut,*(InvM,VecPlusKminBufInc))
    Tbg[2,1]=-dot(Snkminus,*(InvA,Snk))+dot(VecPlusKminBufOut,*(InvM,VecPlusKBufInc))
    Tbg[2,2]=-dot(Snkminus,*(InvA,Snkminus))+dot(VecPlusKminBufOut,*(InvM,VecPlusKminBufInc))

    Sbg=+(Identity,-2*pi*im*Tbg)

   Couplings, Sbg

end

######Input Output relationships**************
function PPvalueRot(HSBarray,freqs,freqidx,l,lprim)

    #TODO, do we have to consider -k as well??
    pts=size(freqs,1)

    integral=convert(ComplexF64,0.0)

    for i in 1:(freqidx-1)
        #integral+=HSBarray[1,l,i]*conj(HSBarray[1,lprim,i])/(freqs[freqidx]-freqs[i])
        #integral+=HSBarray[2,l,i]*conj(HSBarray[2,lprim,i])/(freqs[freqidx]-freqs[i])

        integral+=freqs[i]*HSBarray[1,l,i]*conj(HSBarray[1,lprim,i])/(freqs[freqidx]^2/2-freqs[i]^2/2)
        integral+=freqs[i]*HSBarray[2,l,i]*conj(HSBarray[2,lprim,i])/(freqs[freqidx]^2/2-freqs[i]^2/2)
        #integral=integral/(freqs[freqidx]-freqs[i])

    end

    for i in (freqidx+1):pts
        #integral+=HSBarray[1,l,i]*conj(HSBarray[1,lprim,i])/(freqs[freqidx]-freqs[i])
        #integral+=HSBarray[2,l,i]*conj(HSBarray[2,lprim,i])/(freqs[freqidx]-freqs[i])


        integral+=freqs[i]*HSBarray[1,l,i]*conj(HSBarray[1,lprim,i])/(freqs[freqidx]^2/2-freqs[i]^2/2)
        integral+=freqs[i]*HSBarray[2,l,i]*conj(HSBarray[2,lprim,i])/(freqs[freqidx]^2/2-freqs[i]^2/2)
        #integral=integral/(freqs[freqidx]-freqs[i])
    end
    integral=-integral*abs(freqs[1]-freqs[2])

    #integral*0
end

#Function that builds the 'relaxation matrix'
#dim is the number of system modes we consider in the calculation
#freqsys is an array of size dim that contains the frequency of the system modes
#HSBarray is a 2 \times dim \times N array that contains N system-bath couplings for the range of freq of interest
#HSBarray[0,x,:] corresponds to the couplings for incoming modes from the left, to system mode x+1
#HSBarray[1,y,:] corresponds to the couplings for incoming modes from the right, to system mode y+1

#freqs is the corresponding frequency array
#freqidx is the index within 'freqs' that identifies the freq at which we compute reflection and transmission

function BuildDllRot(dim,freqsys,HSBarray,freqs,freqidx)

    Dll=zeros(ComplexF64,(dim,dim))

    for i in 1:dim

        for j in 1:dim

            Dll[i,j]=PPvalueRot(HSBarray,freqs,freqidx,i,j)
            Dll[i,j]=Dll[i,j]+pi*im*HSBarray[1,i,freqidx]*conj(HSBarray[1,j,freqidx])
            Dll[i,j]=Dll[i,j]+pi*im*HSBarray[2,i,freqidx]*conj(HSBarray[2,j,freqidx])
        end
    end


    for i in 1:dim
        #Dll[i,i]=Dll[i,i]+freqs[freqidx]*(freqs[freqidx]-freqsys[i])
        Dll[i,i]=Dll[i,i]+0.5*(freqs[freqidx]^2-freqsys[i]^2)
    end

    Dll
end

#function to build the IO scattering matrix for a given freq, with index freqidx within array 'freqs'
function BuildSRot(dim,freqsys,HSBarray,freqs,freqidx)

    #Note that this function computes the I/O scattering matrix within the rotating wave-approximation
    #Lendtrot makes a correspondence to the Schrodinger equation to come up with the calculations of
    #the I/O scattering matrix consistent with this approximation. HSBarray must be multiplies
    #by sqrt(\omega_\lambda) before this process as a a result of this scheme!



    DMat=BuildDllRot(dim,freqsys,HSBarray,freqs,freqidx)
    column=zeros(ComplexF64,dim)
    row=zeros(ComplexF64,dim)

    S=zeros(ComplexF64,(2,2))

    #reorganization of data:
    Iden=zeros(ComplexF64,(2,2))
    Iden[1,1]=1
    Iden[2,2]=1


    for m in 1:2

        for mprim in 1:2

            for i in 1:dim
                column[i]=HSBarray[mprim,i,freqidx]
                row[i]=HSBarray[m,i,freqidx]
            end

        S[m,mprim]= dot(row,*(inv(DMat),column))
        S[m,mprim]=-2*pi*im*S[m,mprim]
        end

    end
    S=S+Iden

end





###********************************* DEPRECATED CODE, KEPT FOR THE RECORD*********
#calculation of system-bath couplings, and background scattering matrix, in only one function.
#grid is the array containing the positions over which the dielectric permittivity profile is sampled
#rfs is the array containing the position in space where the dielectric is different than 1
#V is the potential sampled at the regions where the permittivity is different than 1. IMPORTANTLY
#it is assumed that it the potential is already weighted by the extra \epsilon factor that results from the inner
#product and any other relevant scaling factors

#epsilon is the array containing the permittivities that correspond to array "grid"
#Systmods is an array of size gridpts x dim, where gridpts=size(grid) and dim is the number of system modes
#assumed in our calculation
#ModFreqs is an array of size dim that contains the frequencies of the system modes under consideration
#k is the wavenumber at which we compute the system-bath coupling

function calcSysBathPos(grid,rfs,V,epsilon,SystMods,ModFreqs,k)

    kabs=abs(k)
    gridpts,dim=size(SystMods)
    NsystDim=size(rfs,1)


    freeK=(1/sqrt(2*pi*kabs))*exp.(im*kabs*grid) #Using energy-normalized free-partcle states!
    freeminusK=(1/sqrt(2*pi*kabs))*exp.(-im*kabs*grid)
    #Building of all relevant vectors and matrices. Notice this is done such that we avoid redundant
    #calculations as much as possible:

    #vectors...
    Snk=zeros(ComplexF64,dim)
    Snkminus=zeros(ComplexF64,dim)

    VNsysKPlus=zeros(ComplexF64,NsystDim)
    VNsysKMinus=zeros(ComplexF64,NsystDim)



    #Matrices...
    #BMat=zeros(ComplexF64,(dim,dim))
    CMat=zeros(ComplexF64,(NsystDim,dim))
    CMatOut=zeros(ComplexF64,(NsystDim,dim)) #For the calculation of background scatetring matrix
    DMat=zeros(ComplexF64,(dim,NsystDim))

    Tbg=zeros(ComplexF64,(2,2))
    Identity=zeros(ComplexF64,(2,2))
    Identity[1,1]=1.0
    Identity[2,2]=1.0

    #To store final results
    Couplings=zeros(ComplexF64,(2,dim))


    for i in 1:dim
        Snk[i]=MaxwellInner(grid,epsilon,SystMods[:,i],freeK)

        Snkminus[i]=MaxwellInner(grid,epsilon,SystMods[:,i],freeminusK)

    end


    for i in 1:NsystDim

        VNsysKPlus[i]=V[i]*(1/sqrt(2*pi*kabs))*exp(im*kabs*rfs[i])
        VNsysKMinus[i]=V[i]*(1/sqrt(2*pi*kabs))*exp(-im*kabs*rfs[i])
    end

    for i in 1:dim

        dumG0Chi=ActRetGreensPos(grid,SystMods[:,i],rfs,k)
        dumG0Chiconj=conj(dumG0Chi)


        for l in 1:NsystDim

            CMat[l,i]=V[l]*dumG0Chi[l]
            CMatOut[l,i]=V[l]*dumG0Chiconj[l]
        end

    end


    for i in 1:dim

        dumG0Chi=ActRetGreensPos(grid,SystMods[:,i],rfs,k)

        for j in 1:NsystDim

            DMat[i,j]=dumG0Chi[j]*V[j]

        end

    end

    #calculation of InvA and InvM
    InvA=InvGreens(grid,epsilon,SystMods,k)
    InvAOut=conj(InvA)

    InvM=NonOverInvM(grid,epsilon,SystMods,rfs,InvA,V,k)


    #End of vector and matrix calculations *****

    #auxiliary vector to avoid redundant calculations
    VecPlusKBufInc=VNsysKPlus-*(CMat,*(InvA,Snk))            #useful for the calculation of the Tbg matrix
    VecPlusKminBufInc=VNsysKMinus-*(CMat,*(InvA,Snkminus))

    #For computation of background scattering matrix
    VecPlusKBufOut=VNsysKPlus-*(CMatOut,*(InvAOut,Snk))            #useful for the calculation of the Tbg matrix
    VecPlusKminBufOut=VNsysKMinus-*(CMatOut,*(InvAOut,Snkminus))

    VecPlusK=*(InvM,VecPlusKBufInc)
    VecPlusKmin=*(InvM,VecPlusKminBufInc)

    dumKin=*(InvA,*(DMat,VecPlusK))
    dumKinmin=*(InvA,*(DMat,VecPlusKmin))
    Kinetic=*(InvA,Snk)+dumKin
    KineticR=*(InvA,Snkminus)+dumKinmin



    #final step for calculation of coupling... DEPENDING ON THE REPRESENTATION OF BATH COORDINATES, WE WILL
                                                #HAVE A DIFFERENT SCALAR PREFACTOR

    Couplings[1,:]=Kinetic
    Couplings[2,:]=KineticR

    for i in 1:dim
        Couplings[1,i]=0.5*Couplings[1,i]#/sqrt(abs(k)*ModFreqs[i])
        Couplings[2,i]=0.5*Couplings[2,i]#/sqrt(abs(k)*ModFreqs[i])


    end

    #calculation of Tbg matrix, we use the ordering [[++,+-],[-+,--]]. Remember that dot function
    #takes complex conjugate on left vector
    Tbg[1,1]=-dot(Snk,*(InvA,Snk))+dot(VecPlusKBufOut,*(InvM,VecPlusKBufInc))
    Tbg[1,2]=-dot(Snk,*(InvA,Snkminus))+dot(VecPlusKBufOut,*(InvM,VecPlusKminBufInc))
    Tbg[2,1]=-dot(Snkminus,*(InvA,Snk))+dot(VecPlusKminBufOut,*(InvM,VecPlusKBufInc))
    Tbg[2,2]=-dot(Snkminus,*(InvA,Snkminus))+dot(VecPlusKminBufOut,*(InvM,VecPlusKminBufInc))

    Sbg=+(Identity,-2*pi*im*Tbg)

   Couplings, Sbg

end
##*****************SECTION OF FUNCTIONS TO COMPUTE THE INPUT-OUTPUT CONTRIBUTIONS TO SCATTERING******

#HSBarray is a 2 \times dim \times N array that contains N system-bath couplings for the range of freq of interest
#HSBarray[0,x,:] corresponds to the couplings for incoming modes from the left, to system mode x+1
#HSBarray[1,y,:] corresponds to the couplings for incoming modes from the right, to system mode y+1

#freqs is the corresponding frequency array
#freqidx is the index within 'freqs' that identifies the freq at which we compute reflection and transmission

function PPvalue(HSBarray,freqs,freqidx,l,lprim)

    #TODO, do we have to consider -k as well??
    pts=size(freqs,1)

    integral=convert(ComplexF64,0.0)

    for i in 1:(freqidx-1)
        integral+=HSBarray[1,l,i]*conj(HSBarray[1,lprim,i])/(freqs[freqidx]-freqs[i])
        integral+=HSBarray[2,l,i]*conj(HSBarray[2,lprim,i])/(freqs[freqidx]-freqs[i])
        #integral=integral/(freqs[freqidx]-freqs[i])

    end

    for i in (freqidx+1):pts
        integral+=HSBarray[1,l,i]*conj(HSBarray[1,lprim,i])/(freqs[freqidx]-freqs[i])
        integral+=HSBarray[2,l,i]*conj(HSBarray[2,lprim,i])/(freqs[freqidx]-freqs[i])
        #integral=integral/(freqs[freqidx]-freqs[i])
    end
    integral=-integral*abs(freqs[1]-freqs[2])

end

#Function that builds the 'relaxation matrix'
#dim is the number of system modes we consider in the calculation
#freqsys is an array of size dim that contains the frequency of the system modes
#HSBarray is a 2 \times dim \times N array that contains N system-bath couplings for the range of freq of interest
#HSBarray[0,x,:] corresponds to the couplings for incoming modes from the left, to system mode x+1
#HSBarray[1,y,:] corresponds to the couplings for incoming modes from the right, to system mode y+1

#freqs is the corresponding frequency array
#freqidx is the index within 'freqs' that identifies the freq at which we compute reflection and transmission

function BuildDll(dim,freqsys,HSBarray,freqs,freqidx)
    #TODO, do we have to consider -k as well??

    Dll=zeros(ComplexF64,(dim,dim))

    for i in 1:dim

        for j in 1:dim

            Dll[i,j]=PPvalue(HSBarray,freqs,freqidx,i,j)
            Dll[i,j]=Dll[i,j]+pi*im*HSBarray[1,i,freqidx]*conj(HSBarray[1,j,freqidx])
            Dll[i,j]=Dll[i,j]+pi*im*HSBarray[2,i,freqidx]*conj(HSBarray[2,j,freqidx])
        end
    end


    for i in 1:dim
        Dll[i,i]=Dll[i,i]+freqs[freqidx]-freqsys[i]
    end

    Dll
end

#function to build the IO scattering matrix for a given freq, with index freqidx within array 'freqs'
function BuildS(dim,freqsys,HSBarray,freqs,freqidx)

    DMat=BuildDll(dim,freqsys,HSBarray,freqs,freqidx)
    column=zeros(ComplexF64,dim)
    row=zeros(ComplexF64,dim)

    S=zeros(ComplexF64,(2,2))

    #reorganization of data:
    Iden=zeros(ComplexF64,(2,2))
    Iden[1,1]=1
    Iden[2,2]=1


    for m in 1:2

        for mprim in 1:2

            for i in 1:dim
                column[i]=HSBarray[mprim,i,freqidx]
                row[i]=HSBarray[m,i,freqidx]
            end

        S[m,mprim]= dot(row,*(inv(DMat),column))
        S[m,mprim]=-2*pi*im*S[m,mprim]
        end

    end
    S=S+Iden

end

#This function returns a) the range of frequencies (in eVs) spanned by the range defined by the array V (see below)
#b) the 2 x 2 x M array that contains the total scattering matrix, for the range of freqs defined in a)
#and c) the absorption spectrum in the same range of frequencies
#
#INPUT parameters:
#1) Freqs is the array of frequencies considered in the compiutation of system-bath couplings
#2) the path to file that contains the the expectation values of micrcovity amplitude and bath. NOTICE: it is assumed
#an specific order in the information stored, check an example.
#3) two-element array V, such that [V[1],V[2]] denotes the range of INDEXES of the frequency array used for the
#calculation of system-bath couplings.
#4) HSB, the array that contains all the system-bath couplings
#5) ScatMats is the background scattering matrix that needs to be computed before
#6) AbsPath is a N x 2 array, where the (i,1)th entry corresponds to path to the time-dep. coherence \langle i | O_{exc}(t) |ref\rangle
# sampling points; where O_exc(t) is the total polarization operator and | i\rangle is an auxiliary state.
#Similarly, the (i,2)th entry is the same but substituting O_{exc}(t) -> q_{ph} is the photonic coordinate.
#7) mu, is the single light-matter coupling
#8) we, is the electronic transition energy

function LinSpect(Freqs,Path,V,HSB,ScatMats,AbsPath,mu,we)
    #useful parameters:
    CovToEv=0.197*2
    
    
    #****FIRST PHASE: read the expectation values of cavity and bath mode and Fourier transform it****
    AmpDum=readdlm(Path);
    

    Timepts,columns=size(AmpDum)

    
    MicAmpMol=zeros(ComplexF64,(Timepts,2))
    
    BathAmpMol=zeros(ComplexF64,(Timepts,2))

    

    for i in 1:Timepts
        scale=AmpDum[i,2]

        MicAmpMol[i,1]=AmpDum[i,1] #Time coordinate
        MicAmpMol[i,2]=AmpDum[i,3]*scale #real part of coordinate
        MicAmpMol[i,2]+=im*AmpDum[i,4]*scale


        BathAmpMol[i,1]=AmpDum[i,1] #Time coordinate
        BathAmpMol[i,2]=AmpDum[i,5]*scale #real part of coordinate
        BathAmpMol[i,2]+=im*AmpDum[i,6]*scale

        
    end




    deltaT=abs(MicAmpMol[1,1]-MicAmpMol[2,1]);
    Ttime=abs(MicAmpMol[1,1]-MicAmpMol[Timepts,1]) #This is the total time of the simulation, in fs!!
    
    FouMicMol=ifft(MicAmpMol[:,2])*length(MicAmpMol[:,2])*deltaT

    FreqsFFTMol=0.658*fftfreq(length(MicAmpMol[:,2]))*2*pi/deltaT ;#Here we take advantage of the equal size along all dimensions

    FreqsFFTMol=fftshift(FreqsFFTMol)
    FouMicMol=fftshift(FouMicMol);
    
    FouBathMol=ifft(BathAmpMol[:,2])*length(MicAmpMol[:,2])*deltaT
    FouBathMol=fftshift(FouBathMol);
    
    #interpolation...
    spRMicMol=Spline1D(real(FreqsFFTMol), real(FouMicMol); w=ones(length(FreqsFFTMol)), k=3, bc="nearest", s=0.0);
    spImMicMol=Spline1D(real(FreqsFFTMol),imag(FouMicMol); w=ones(length(FreqsFFTMol)), k=3, bc="nearest", s=0.0);
    spRBathMol=Spline1D(real(FreqsFFTMol), real(FouBathMol); w=ones(length(FreqsFFTMol)), k=3, bc="nearest", s=0.0);
    spImBathMol=Spline1D(real(FreqsFFTMol), imag(FouBathMol); w=ones(length(FreqsFFTMol)), k=3, bc="nearest", s=0.0);
    
        
    #lets consider a range of frequencies for which we know the system-bath coupling constants
    FRangeIdx=V[1]:1:V[2]

    Rangepts=length(FRangeIdx)
    FRangeIO=Freqs[V[1]:V[2]].^2 *CovToEv/2 #Freqs in eVs!!! These parameters are already defined in the molecular-free
    #section

    #calculation "by hand" of the IO scattering matrix for each frequency
    SIO_NumMol=zeros(ComplexF64,(2,2,Rangepts))


    for i in 1:Rangepts
       #retrieve the system-bath coupling for each frequency, as well as the relative phase between the degenerate modes
        #conversion to sqrt(eV), those are the units of this guy...
        CoupLeft=sqrt(CovToEv)*HSB[1,1,FRangeIdx[i]] #left
        CoupRight=sqrt(CovToEv)*HSB[2,1,FRangeIdx[i]]; #right

        PsiLeft=angle(CoupLeft)
        PsiRight=angle(CoupRight)

        ELR=exp(-im*(PsiLeft-PsiRight)) #difference in phases of couplings
        AuxRatio=abs(CoupRight)/abs(CoupLeft)
        #ratio of bath and microcavity amplitudes...We assume that the bath that drives the system is always from the left
        #such that the we can scale that result appropiately to get the corresponding result for the right bath
        Ratio=(evaluate(spRMicMol,FRangeIO[i])+im*evaluate(spImMicMol,FRangeIO[i]))/(evaluate(spRBathMol,FRangeIO[i])+im*evaluate(spImBathMol,FRangeIO[i]))

        SIO_NumMol[1,1,i]=1-2*pi*im*abs(CoupLeft)*Ratio
        SIO_NumMol[1,2,i]=-2*pi*im*abs(CoupLeft)*AuxRatio*ELR*Ratio
        SIO_NumMol[2,1,i]=-2*pi*im*abs(CoupRight)*conj(ELR)*Ratio
        SIO_NumMol[2,2,i]=1-2*pi*im*abs(CoupRight)*AuxRatio*Ratio



    end
    
    STot_NumMol=zeros(ComplexF64,(2,2,Rangepts))

    for i in 1:Rangepts

        STot_NumMol[:,:,i]=*(ScatMats[:,:,FRangeIdx[i]],SIO_NumMol[:,:,i])

    end
    
    #******SECOND PHASE: Read the coherences introduced as input to compute absorption******
    #extract dimensions of path array:
    NAuxStates,=size(AbsPath)
    
    init=1
    
    endT=Timepts

    WRange=10 #Frequency range goes from -WRange/2 to WRange/2, in eV

    SampTime=2*pi*0.658/WRange

    #sampling rate in fs
    SampRate=Int64(div(2*pi*0.658/WRange,deltaT,RoundDown))

    #total number of points:
    TotPts=endT-init+1

    #number of points to sample:
    SampPts=Int64(div(TotPts,SampRate,RoundDown))
    
    FreqsAbs=0.658*fftfreq(SampPts)*2*pi/SampTime
    FreqsAbs=fftshift(FreqsAbs)
    
    TotAbs=zeros(Float64,SampPts)
       
    for i in 1:NAuxStates
        
        PolCoh=readdlm(AbsPath[i,1])
        PhotCoh=readdlm(AbsPath[i,2])

        FilteredPhot=zeros(ComplexF64,SampPts)
        FilteredExc=zeros(ComplexF64,SampPts)

        for j in 1:SampPts
            
            FilteredPhot[j]=conj(PhotCoh[init+j*SampRate-1,2]+im*PhotCoh[init+j*SampRate-1,3])
            FilteredExc[j]=PolCoh[init+j*SampRate-1,2]+im*PolCoh[init+j*SampRate-1,3]

        end

        FouAbsPhot=ifft(FilteredPhot)*TotPts*deltaT/0.658
        FouAbsExc=ifft(FilteredExc)*TotPts*deltaT/0.658
        
        FouAbsPhot=fftshift(FouAbsPhot)
        FouAbsExc=fftshift(FouAbsExc);


        Product=conj(FouAbsPhot).*FouAbsExc
        
        TotAbs+=-imag(Product)*2*sqrt(2)*mu*we


        
        
    end
    
    #*****THIRD PHASE: finally, we interpolate the bath-mode data and calculate the absorption****
    
    # we estimate the indexes of the array that correspond to the frequencies contained in FRangeIO
    deltaFreqBath=abs(FreqsFFTMol[2]-FreqsFFTMol[1])
    initFreqBath=Int64(div(Timepts,2,RoundDown))+ Int64(div(FRangeIO[1],deltaFreqBath,RoundDown))
    finFreqBath=Int64(div(Timepts,2,RoundDown))+ Int64(div(FRangeIO[length(FRangeIO)],deltaFreqBath,RoundDown))
    
    
    spInput=Spline1D(FreqsFFTMol[initFreqBath:finFreqBath], abs.(FouBathMol[initFreqBath:finFreqBath]).^2; w=ones(length(FreqsFFTMol[initFreqBath:finFreqBath])),
    k=3, bc="nearest", s=0.0);
    
    #interpolation of the absorption within the range of frequencies of FRangeIO
    deltaFreqAbs=abs(FreqsAbs[2]-FreqsAbs[1])
    initFreqAbs=Int64(div(SampPts,2,RoundDown))+ Int64(div(FRangeIO[1],deltaFreqAbs,RoundDown))
    finFreqAbs=Int64(div(SampPts,2,RoundDown))+ Int64(div(FRangeIO[length(FRangeIO)],deltaFreqAbs,RoundDown))
    

    NormAbs=zeros(Float64,length(FreqsAbs[initFreqAbs:finFreqAbs]))
    
    RenormAbs=zeros(Float64,length(FreqsAbs[initFreqAbs:finFreqAbs]))

    for i in 1:length(FreqsAbs[initFreqAbs:finFreqAbs])

        NormAbs[i]=TotAbs[initFreqAbs+i-1]
        RenormAbs[i]=(FreqsAbs[initFreqAbs+i-1]*(evaluate(spInput,FreqsAbs[initFreqAbs+i-1]))/0.658)
        #NormAbs[i]=TotalAbs[150+i-1]/abs(evaluate(spInput,FreqsTest[150+i-1]))

    end
    
    #Finally, we compute the absorption spectrum using the fine-grained frequency grid of FRangeIO
    spAbs=Spline1D(FreqsAbs[initFreqAbs:finFreqAbs], NormAbs; w=ones(length(FreqsAbs[initFreqAbs:finFreqAbs])), k=3, bc="nearest", s=0.0)
    spRenormAbs=Spline1D(FreqsAbs[initFreqAbs:finFreqAbs], RenormAbs; w=ones(length(FreqsAbs[initFreqAbs:finFreqAbs])), k=3, bc="nearest", s=0.0)

    AdjustedAbs=zeros(Float64,length(FRangeIO))
    AdjustedRenormAbs=zeros(Float64,length(FRangeIO))

    for i in 1:length(FRangeIO)

        AdjustedAbs[i]= evaluate(spAbs,FRangeIO[i])
        AdjustedRenormAbs[i]=evaluate(spRenormAbs,FRangeIO[i])
        

    end
    
    #Based on our assumptions for the simple model, we can compute the ratio of the reflection amplitudes without molecule(s)
    # and that with molecules:
    
    #rho=
    
    
    FRangeIO,STot_NumMol,AdjustedAbs,AdjustedRenormAbs
    #FreqsAbs,TotAbs

    
    
    
    
end

#I am defining this function to compute the linear spectra using the newest approach that does not
#rely on the external photon mode
#INPUT parameters:
#1) Freqs is the array of frequencies considered in the compiutation of system-bath couplings
#2) Init_Amp refers to the initial (time=0) expectation value of the effective amplitude of the external driving
#mode, which corresponds to the microcavity amplitude at t=0 divided by the absolute value of
#of the matrix tunnelling element at the effetive frequency of the mode :S
#3) the path to file that contains the the expectation values of micrcovity amplitude and bath. NOTICE: it is assumed
#an specific order in the information stored, check an example.
#4) two-element array V, such that [V[1],V[2]] denotes the range of INDEXES of the frequency array used for the
#calculation of system-bath couplings.
#5) HSB, the array that contains all the system-bath couplings
#6) ScatMats is the background scattering matrix that needs to be computed before
#7) AbsPath is a N x 2 array, where the (i,1)th entry corresponds to path to the time-dep. coherence \langle i | O_{exc}(t) |ref\rangle
# sampling points; where O_exc(t) is the total polarization operator and | i\rangle is an auxiliary state.
#Similarly, the (i,2)th entry is the same but substituting O_{exc}(t) -> q_{ph} is the photonic coordinate.
#8) mu, is the single light-matter coupling
#9) we, is the electronic transition energy



function LinSpectNew(Freqs,Init_Amp,Path,V,HSB,ScatMats,AbsPath,mu,we)
    #useful parameters:
    CovToEv=0.197*2
    
    
    
    #****FIRST PHASE: read the expectation values of cavity and bath mode and Fourier transform it****
    AmpDum=readdlm(Path);
    

    Timepts,columns=size(AmpDum)

    
    MicAmpMol=zeros(ComplexF64,(Timepts,2))
    
    #BathAmpMol=zeros(ComplexF64,(Timepts,2))

    

    for i in 1:Timepts
        scale=AmpDum[i,2]

        MicAmpMol[i,1]=AmpDum[i,1] #Time coordinate
        MicAmpMol[i,2]=AmpDum[i,3]*scale #real part of coordinate
        MicAmpMol[i,2]+=im*AmpDum[i,4]*scale


        #BathAmpMol[i,1]=AmpDum[i,1] #Time coordinate
        #BathAmpMol[i,2]=AmpDum[i,5]*scale #real part of coordinate
        #BathAmpMol[i,2]+=im*AmpDum[i,6]*scale

        
    end




    deltaT=abs(MicAmpMol[1,1]-MicAmpMol[2,1]);
    Ttime=abs(MicAmpMol[1,1]-MicAmpMol[Timepts,1]) #This is the total time of the simulation, in fs!!
    
    FouMicMol=ifft(MicAmpMol[:,2])*length(MicAmpMol[:,2])*deltaT

    FreqsFFTMol=0.658*fftfreq(length(MicAmpMol[:,2]))*2*pi/deltaT ;#Here we take advantage of the equal size along all dimensions

    FreqsFFTMol=fftshift(FreqsFFTMol)
    FouMicMol=fftshift(FouMicMol);
    
    #FouBathMol=ifft(BathAmpMol[:,2])*length(MicAmpMol[:,2])*deltaT
    #FouBathMol=fftshift(FouBathMol);
    
    #interpolation...
    spRMicMol=Spline1D(real(FreqsFFTMol), real(FouMicMol); w=ones(length(FreqsFFTMol)), k=3, bc="nearest", s=0.0);
    spImMicMol=Spline1D(real(FreqsFFTMol),imag(FouMicMol); w=ones(length(FreqsFFTMol)), k=3, bc="nearest", s=0.0);
    #spRBathMol=Spline1D(real(FreqsFFTMol), real(FouBathMol); w=ones(length(FreqsFFTMol)), k=3, bc="nearest", s=0.0);
    #spImBathMol=Spline1D(real(FreqsFFTMol), imag(FouBathMol); w=ones(length(FreqsFFTMol)), k=3, bc="nearest", s=0.0);
    
        
    #lets consider a range of frequencies for which we know the system-bath coupling constants
    FRangeIdx=V[1]:1:V[2]

    Rangepts=length(FRangeIdx)
    FRangeIO=Freqs[V[1]:V[2]].^2 *CovToEv/2 #Freqs in eVs!!! These parameters are already defined in the molecular-free
    #section

    #calculation "by hand" of the IO scattering matrix for each frequency
    SIO_NumMol=zeros(ComplexF64,(2,2,Rangepts))


    for i in 1:Rangepts
       #retrieve the system-bath coupling for each frequency, as well as the relative phase between the degenerate modes
        #conversion to sqrt(eV), those are the units of this guy...
        CoupLeft=sqrt(CovToEv)*HSB[1,1,FRangeIdx[i]] #left
        CoupRight=sqrt(CovToEv)*HSB[2,1,FRangeIdx[i]]; #right

        PsiLeft=angle(CoupLeft)
        PsiRight=angle(CoupRight)

        ELR=exp(-im*(PsiLeft-PsiRight)) #difference in phases of couplings
        AuxRatio=abs(CoupRight)/abs(CoupLeft)
        #ratio of bath and microcavity amplitudes...We assume that the bath that drives the system is always from the left
        #such that the we can scale that result appropiately to get the corresponding result for the right bath
        #Ratio=(evaluate(spRMicMol,FRangeIO[i])+im*evaluate(spImMicMol,FRangeIO[i]))/(evaluate(spRBathMol,FRangeIO[i])+im*evaluate(spImBathMol,FRangeIO[i]))
        Ratio=(evaluate(spRMicMol,FRangeIO[i])+im*evaluate(spImMicMol,FRangeIO[i]))/Init_Amp


        #SIO_NumMol[1,1,i]=1-2*pi*im*abs(CoupLeft)*Ratio
        #SIO_NumMol[1,2,i]=-2*pi*im*abs(CoupLeft)*AuxRatio*ELR*Ratio
        #SIO_NumMol[2,1,i]=-2*pi*im*abs(CoupRight)*conj(ELR)*Ratio
        #SIO_NumMol[2,2,i]=1-2*pi*im*abs(CoupRight)*AuxRatio*Ratio
        
        SIO_NumMol[1,1,i]=1-2*pi*abs(CoupLeft)*Ratio
        SIO_NumMol[1,2,i]=-2*pi*abs(CoupLeft)*AuxRatio*ELR*Ratio
        SIO_NumMol[2,1,i]=-2*pi*abs(CoupRight)*conj(ELR)*Ratio
        SIO_NumMol[2,2,i]=1-2*pi*abs(CoupRight)*AuxRatio*Ratio



    end
    
    STot_NumMol=zeros(ComplexF64,(2,2,Rangepts))

    for i in 1:Rangepts

        STot_NumMol[:,:,i]=*(ScatMats[:,:,FRangeIdx[i]],SIO_NumMol[:,:,i])

    end
    
    #******SECOND PHASE: Read the coherences introduced as input to compute absorption******
    #extract dimensions of path array:
    NAuxStates,=size(AbsPath)
    
    init=1
    
    endT=Timepts

    WRange=10 #Frequency range goes from -WRange/2 to WRange/2, in eV

    SampTime=2*pi*0.658/WRange

    #sampling rate in fs
    SampRate=Int64(div(2*pi*0.658/WRange,deltaT,RoundDown))

    #total number of points:
    TotPts=endT-init+1

    #number of points to sample:
    SampPts=Int64(div(TotPts,SampRate,RoundDown))
    
    FreqsAbs=0.658*fftfreq(SampPts)*2*pi/SampTime
    FreqsAbs=fftshift(FreqsAbs)
    
    TotAbs=zeros(Float64,SampPts)
    
    #Array to store the different Fourier transforms of coherences
    FouAbsPhotMaster=zeros(ComplexF64,(NAuxStates,SampPts))
    FouAbsExcMaster=zeros(ComplexF64,(NAuxStates,SampPts))
    
    for i in 1:NAuxStates
        
        PolCoh=readdlm(AbsPath[i,1])
        PhotCoh=readdlm(AbsPath[i,2])

        FilteredPhot=zeros(ComplexF64,SampPts)
        FilteredExc=zeros(ComplexF64,SampPts)

        for j in 1:SampPts
            
            FilteredPhot[j]=conj(PhotCoh[init+j*SampRate-1,2]+im*PhotCoh[init+j*SampRate-1,3])
            FilteredExc[j]=PolCoh[init+j*SampRate-1,2]+im*PolCoh[init+j*SampRate-1,3]

        end

        FouAbsPhot=ifft(FilteredPhot)*TotPts*deltaT/0.658
        FouAbsExc=ifft(FilteredExc)*TotPts*deltaT/0.658
        
        FouAbsPhot=fftshift(FouAbsPhot)
        FouAbsExc=fftshift(FouAbsExc);
        
        FouAbsPhotMaster[i,:]=FouAbsPhot
        FouAbsExcMaster[i,:]=FouAbsExc


        Product=conj(FouAbsPhot).*FouAbsExc
        
        TotAbs+=-imag(Product)*2*sqrt(2)*mu*we


        
        
    end
    
    #*****THIRD PHASE: finally, we interpolate the bath-mode data and calculate the absorption****
    
    # we estimate the indexes of the array that correspond to the frequencies contained in FRangeIO
    deltaFreqBath=abs(FreqsFFTMol[2]-FreqsFFTMol[1])
    
    #initFreqBath=Int64(div(Timepts,2,RoundDown))+ Int64(div(FRangeIO[1],deltaFreqBath,RoundDown))
    #finFreqBath=Int64(div(Timepts,2,RoundDown))+ Int64(div(FRangeIO[length(FRangeIO)],deltaFreqBath,RoundDown))
    
    
    #spInput=Spline1D(FreqsFFTMol[initFreqBath:finFreqBath], abs.(FouBathMol[initFreqBath:finFreqBath]).^2; w=ones(length(FreqsFFTMol[initFreqBath:finFreqBath])),
    #k=3, bc="nearest", s=0.0);
    
    #interpolation of the absorption within the range of frequencies of FRangeIO
    deltaFreqAbs=abs(FreqsAbs[2]-FreqsAbs[1])
    initFreqAbs=Int64(div(SampPts,2,RoundDown))+ Int64(div(FRangeIO[1],deltaFreqAbs,RoundDown))
    finFreqAbs=Int64(div(SampPts,2,RoundDown))+ Int64(div(FRangeIO[length(FRangeIO)],deltaFreqAbs,RoundDown))
    

    NormAbs=zeros(Float64,length(FreqsAbs[initFreqAbs:finFreqAbs]))
    
    #RenormAbs=zeros(Float64,length(FreqsAbs[initFreqAbs:finFreqAbs]))

    for i in 1:length(FreqsAbs[initFreqAbs:finFreqAbs])

        NormAbs[i]=TotAbs[initFreqAbs+i-1]
        #RenormAbs[i]=(FreqsAbs[initFreqAbs+i-1]*(evaluate(spInput,FreqsAbs[initFreqAbs+i-1]))/0.658)
        #NormAbs[i]=TotalAbs[150+i-1]/abs(evaluate(spInput,FreqsTest[150+i-1]))

    end
    
    #Finally, we compute the absorption spectrum using the fine-grained frequency grid of FRangeIO
    spAbs=Spline1D(FreqsAbs[initFreqAbs:finFreqAbs], NormAbs; w=ones(length(FreqsAbs[initFreqAbs:finFreqAbs])), k=3, bc="nearest", s=0.0)
    #spRenormAbs=Spline1D(FreqsAbs[initFreqAbs:finFreqAbs], RenormAbs; w=ones(length(FreqsAbs[initFreqAbs:finFreqAbs])), k=3, bc="nearest", s=0.0)

    AdjustedAbs=zeros(Float64,length(FRangeIO))
    #AdjustedRenormAbs=zeros(Float64,length(FRangeIO))

    for i in 1:length(FRangeIO)

        AdjustedAbs[i]= evaluate(spAbs,FRangeIO[i])
        #AdjustedRenormAbs[i]=evaluate(spRenormAbs,FRangeIO[i])
        

    end
    
    #Based on our assumptions for the simple model, we can compute the ratio of the reflection amplitudes without molecule(s)
    # and that with molecules:
    
    #rho=
    
    
    FRangeIO,STot_NumMol,AdjustedAbs,spRMicMol,spImMicMol, FouAbsPhotMaster, FouAbsExcMaster, FreqsAbs#,AdjustedRenormAbs
    #FreqsAbs,TotAbs

    
    
    
    
end



