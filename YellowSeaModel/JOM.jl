# main program

include("readtopo.jl")
include("eddyVertical.jl")
include("eddyHorizontal.jl")
include("advection.jl")
#include("temperature.jl")
#include("salinity.jl")
#include("state.jl")
#include("age.jl")
include("hydrostatic.jl")
include("momentadv.jl")
include("bottomShear.jl")
include("surfaceWindShear.jl")
include("momentdif.jl")
include("poisson.jl")
include("writing.jl")

# for parallel computing(not used)
#@everywhere using Distributed
#@everywhere using DistributedArrays
#@everywhere using DistributedArrays.SPMD

# initialize
	nx = 85
	ny = 87
	nz = 26
	dt = 30.0			# time step
	g = 9.81			# gravity acceleration
	rho0 = 1028.0		# background rho
	warmup = 0.0		# warm up parameter
	drag = 0.0015		# bottom friction rate
	slip = 1.0			# slip-nonslip condition indicator
	fbeta = 2.2e-11		# coriolis parameter
	omega = 1.4			# S.O.R parameter
	epsilon = 0.01		# possion accuracy parameter

	periods = 16.0				# total running time (days)
	outinter = 3600.0			# output interval (hours)
	nper = 1				# iteration default number
	
	# tidal forcing 4 points
	fele = 0.0
	fele1 = 0.0
	fele2 = 0.0
	fele3 = 0.0
	fele4 = 0.0
	dele = 0.0
	
	# initializing variables
	dx = zeros(Float64, nx+2,ny+2)
	dy = zeros(Float64, nx+2,ny+2)
	dz = zeros(Float64, nz+2)
	f = zeros(Float64, ny+2)
	ele = zeros(Float64, nx+2,ny+2)
	depth = zeros(Float64, nx+2,ny+2)
	ib = zeros(Int64, nx+2,ny+2)
	land = zeros(Int64, nx+2,ny+2)
	dry = zeros(Int64, nx+2,ny+2,nz+2)
	
	# Velocities
	u = zeros(Float64, nx+2,ny+2,nz+2)
	v = zeros(Float64, nx+2,ny+2,nz+2)
	w = zeros(Float64, nx+2,ny+2,nz+2)
	un = zeros(Float64, nx+2,ny+2,nz+2)
	vn = zeros(Float64, nx+2,ny+2,nz+2)
	wn = zeros(Float64, nx+2,ny+2,nz+2)
	umm = zeros(Float64, nz+1)
	vmm = zeros(Float64, nz+1)
	
	# Eddy Viscosity Variables
	kz = zeros(Float64, nx+2,ny+2,nz+2)
	az = zeros(Float64, nx+2,ny+2,nz+2)
	kh = zeros(Float64, nx+2,ny+2,nz+2)
	ah = zeros(Float64, nx+2,ny+2,nz+2)
	
	# Advection Variables
	t = zeros(Float64, nx+2,ny+2,nz+2)
	tn = zeros(Float64, nx+2,ny+2,nz+2)
	s = zeros(Float64, nx+2,ny+2,nz+2)
	sn = zeros(Float64, nx+2,ny+2,nz+2)
	rho = zeros(Float64, nx+2,ny+2,nz+2)
	rhon = zeros(Float64, nx+2,ny+2,nz+2)
	age = zeros(Float64, nx+2,ny+2,nz+2)
	agen = zeros(Float64, nx+2,ny+2,nz+2)
	B = zeros(Float64, nx+2,ny+2,nz+2)
	BN = zeros(Float64, nx+2,ny+2,nz+2)
	CuPos = zeros(Float64, nx+2,ny+2,nz+2)
	CuNeg = zeros(Float64, nx+2,ny+2,nz+2)
	CvPos = zeros(Float64, nx+2,ny+2,nz+2)
	CvNeg = zeros(Float64, nx+2,ny+2,nz+2)
	CwPos = zeros(Float64, nx+2,ny+2,nz+2)
	CwNeg = zeros(Float64, nx+2,ny+2,nz+2)
	
	# Pressure Variables
	p = zeros(Float64, nx+2,ny+2,nz+2)
	q = zeros(Float64, nx+2,ny+2,nz+2)
	dq = zeros(Float64, nx+2,ny+2,nz+2)
	
	# Momentun Variables
	advx = zeros(Float64, nx+2,ny+2,nz+2)
	advy = zeros(Float64, nx+2,ny+2,nz+2)
	advz = zeros(Float64, nx+2,ny+2,nz+2)
	diffu = zeros(Float64, nx+2,ny+2,nz+2)
	diffv = zeros(Float64, nx+2,ny+2,nz+2)
	diffw = zeros(Float64, nx+2,ny+2,nz+2)
	
	# Stress Variables
	taubu = zeros(Float64, nx+2,ny+2)
	taubv = zeros(Float64, nx+2,ny+2)
	taux = zeros(Float64, nx+2,ny+2)
	tauy = zeros(Float64, nx+2,ny+2)
	
	# S.O.R. Variables
	du = zeros(Float64, nx+2,ny+2,nz+2)
	dv = zeros(Float64, nx+2,ny+2,nz+2)
	dw = zeros(Float64, nx+2,ny+2,nz+2)
	ustar = zeros(Float64, nx+2,ny+2,nz+2)
	vstar = zeros(Float64, nx+2,ny+2,nz+2)
	wstar = zeros(Float64, nx+2,ny+2,nz+2)
	qstar = zeros(Float64, nx+2,ny+2,nz+2)
	atop = zeros(Float64, nx+2,ny+2,nz+2)
	abot = zeros(Float64, nx+2,ny+2,nz+2)
	ae = zeros(Float64, nx+2,ny+2,nz+2)
	aw = zeros(Float64, nx+2,ny+2,nz+2)
	as = zeros(Float64, nx+2,ny+2,nz+2)
	an = zeros(Float64, nx+2,ny+2,nz+2)
	atotal = zeros(Float64, nx+2,ny+2,nz+2)
	uintegral = zeros(Float64, nx+2,ny+2)
	vintegral = zeros(Float64, nx+2,ny+2)
	drdx = zeros(Float64, nx+2,ny+2)
	drdy = zeros(Float64, nx+2,ny+2)
	drdz = zeros(Float64, nz+2)
	nt = zeros(Int64, 1)	# dtin/dt
	np = zeros(Int64, 1)	# iteration number at break
	dtin = zeros(Int64,1)	# exterior time step
	runtime = zeros(Float64,1)	# total runtime
	rundays = zeros(Float64,1)	# total rundays
	
	# read tidal forcing files
	fname1 = "./input/SouthWest.csv"
	fname2 = "./input/SouthEast.csv"
	fname3 = "./input/EastSouth.csv"
	fname4 = "./input/EastNorth.csv"
	ftide1 = readcsv(fname1)
	ftide2 = readcsv(fname2)
	ftide3 = readcsv(fname3)
	ftide4 = readcsv(fname4)
		
		dtin[1] = dt*1
		nt[1] = trunc(Int64,dtin[1]/dt)
	
		@. dx = 10000.0
		@. dy = 10000.0
		@. dz = 5.0
		@. f = 1.e-4
		#@. t = 4.0
		#@. s = 30.0
		#@. age = 0.0
		
		# coriolis parameter setup
		@simd for j = 1:ny+2
			f[j] = 1.e-4 + fbeta*real(j)*dy[2,j]
		end
		
		# read topography from file
		ReadTopo(depth,ele,ib,dx,dy,land,dry)
		# open rho.csv
		fname = "./output/rho.csv"
		csv1 = open(fname,"w+")
		# open u.csv
		fname = "./output/u.csv"
		csv2 = open(fname,"w+")
		# open v.csv
		fname = "./output/v.csv"
		csv3 = open(fname,"w+")
		# open iterations number n.csv
		fname = "./output/w.csv"
		csv4 = open(fname,"w+")
		# open t.csv
		fname = "./output/t.csv"
		csv5 = open(fname,"w+")
		# open s.csv
		fname = "./output/s.csv"
		csv6 = open(fname,"w+")
		# open age.csv
		fname = "./output/age.csv"
		csv7 = open(fname,"w+")
		# open q.csv
		fname = "./output/q.csv"
		csv8 = open(fname,"w+")
		
		write(csv1, ["i, j, rho_1 ~ rho_k \n"])
		write(csv2, ["i, j, u_1 ~ u_k \n"])
		write(csv3, ["i, j, v_1 ~ v_k \n"])
		write(csv4, ["i, j, w_1 ~ w_k \n"])
		write(csv5, ["i, j, t_1 ~ t_k \n"])
		write(csv6, ["i, j, s_1 ~ s_k \n"])
		write(csv7, ["i, j, age_1 ~ age_k \n"])
		write(csv8, ["i, j, q_1 ~ q_k \n"])
		
		@. drdx = 1.0/(rho0*dx)
		@. drdy = 1.0/(rho0*dy)
		@. drdz = 1.0/(rho0*dz)
	
		@. depth = max(depth,0)
		WetDry(nx,ny,nz,u,v,w,q,depth)
		
		# reset the poisson equation parameters
		@fastmath @inbounds for i = 2:nx+1, j = 2:ny+1, k = 2:ib[i,j]
			atop[i,j,k] = dx[i,j]/dz[k]
			abot[i,j,k] = dx[i,j]/dz[k]
			ae[i,j,k] = dz[k]/dx[i,j]
			aw[i,j,k] = dz[k]/dx[i,j]
			an[i,j,k] = dz[k]*dx[i,j]/(dy[i,j]*dy[i,j])
			as[i,j,k] = dz[k]*dx[i,j]/(dy[i,j]*dy[i,j])
					if (dry[i+1,j,k] != 0) ae[i,j,k] = 0.0 end
					if (i == nx+1) ae[i,j,k] = 0.0 end
					if (dry[i-1,j,k] != 0) aw[i,j,k] = 0.0 end
					if (i == 2) aw[i,j,k] = 0.0 end
					if (dry[i,j+1,k] != 0) an[i,j,k] = 0.0 end
					if (j == ny+1) an[i,j,k] = 0.0 end
					if (dry[i,j-1,k] != 0) as[i,j,k] = 0.0 end
					if (j == 2) as[i,j,k] = 0.0 end
					if (dry[i,j,k+1] != 0) abot[i,j,k] = 0.0 end
					if (dry[i,j,k-1] != 0) atop[i,j,k] = 0.0 end
			atotal[i,j,k] = abot[i,j,k] + atop[i,j,k] + ae[i,j,k] + aw[i,j,k] + an[i,j,k] + as[i,j,k]
		end
		
	# wind-stress forcing
	taux .= 0.0
	tauy .= 0.0
			
	# runtime set
	ntotal = trunc(Int64,periods*24*60*60/dtin[1])
	nout = trunc(Int64,outinter/dtin[1])
		
	# output initial values	
	# write h.csv
		fname = "./output/h.csv"
		csvfile = open(fname,"w")
		write(csvfile, ["i, j, dx, dy, depth, land, ib\n"])
		@inbounds for j = 1:ny+2, i = 1:nx+2
				write(csvfile, [string.(i),",",string.(j),",",string.(dx[i,j]),",",string.(dy[i,j]),",",string.(depth[i,j]),",",string.(land[i,j]),",",string.(ib[i,j]),"\n"])
		end
		close(csvfile)
	# write dry.csv
		fname = "./output/dry.csv"
		csvfile = open(fname,"w")
		write(csvfile, ["i, j, k, depth, ib, land, dry\n"])
		@inbounds for k = 1:nz+2, j = 1:ny+2, i = 1:nx+2
				write(csvfile, [string.(i),",",string.(j),",",string.(k),",",string.(depth[i,j]),",",string.(ib[i,j]),",",string.(land[i,j]),",",string.(dry[i,j,k]),"\n"])
		end
		close(csvfile)
	# write initial rho.csv	
		@inbounds for j = 2:ny+1, i = 2:nx+1
					write(csv1, [string.(i),",",string.(j),",",string.(rho[i,j,2:nz+1],","),"\n"])
					write(csv8, [string.(i),",",string.(j),",",string.(q[i,j,2:nz+1],","),"\n"])
			end
		
		close(csv1)
		close(csv2)
		close(csv3)
		close(csv4)
		close(csv5)
		close(csv6)
		close(csv7)
		close(csv8)
		
		EddyVertical(nx,ny,nz,dz,u,v,rho,kz,az,rho0,g) 
		EddyHorizontal(nx,ny,nz,dx,dy,u,v,w,ah,kh) 
		BottomShear(nx,ny,ib,drag,u,v,taubu,taubv)
		Hydrostatic(nx,ny,nz,dz,g,p,rho)
	
############### simulation loop start
	@inbounds for n = 1:ntotal
				
			runtime[1] = runtime[1] + dtin[1]
			rundays[1] = runtime[1]/86400.0
			
			@. depth[2:nx+1,2:ny+1] = max( -ele[2:nx+1,2:ny+1] + q[2:nx+1,2:ny+1,1]/(g*rho0), 0)
			 WetDry(nx,ny,nz,u,v,w,q,depth)
				
		# flow dynamic calculate
			 EddyVertical(nx,ny,nz,dz,u,v,rho,kz,az,rho0,g) 
			 EddyHorizontal(nx,ny,nz,dx,dy,u,v,w,ah,kh) 
			 BottomShear(nx,ny,ib,drag,u,v,taubu,taubv) 
			 SurfaceWindShear(nx,ny,dz,warmup,rho0,az,u,v,taux,tauy) 
			#
			 #Age(nx,ny,nz,dry,dx,dy,dz,dt,dtin,age,agen,B,BN,u,v,w,ah,kz,ib,CuPos,CuNeg,CvPos,CvNeg,CwPos,CwNeg)
 			 #Temperature(nx,ny,nz,dry,dx,dy,dz,dtin,t,tn,B,BN,u,v,w,ah,kz,CuPos,CuNeg,CvPos,CvNeg,CwPos,CwNeg) 
 			 #Salinity(nx,ny,nz,dry,dx,dy,dz,dtin,s,sn,B,BN,u,v,w,ah,kz,CuPos,CuNeg,CvPos,CvNeg,CwPos,CwNeg) 
 			 #State(t,s,rho,rho0)
 			 Hydrostatic(nx,ny,nz,dz,g,p,rho) 
			 MomentumAdvection(nx,ny,nz,dx,dy,dz,dt,dtin,u,v,w,B,BN,advx,advy,advz,CuPos,CuNeg,CvPos,CvNeg,CwPos,CwNeg) 
			 MomentumDiffusion(nx,ny,nz,dx,dy,dz,dtin,ib,u,v,w,ah,az,dry,land,slip,taubu,diffu,taubv,diffv,diffw) 
			
		# tidal forcing	
			lt=length(ftide1)
			for i = 1:lt-1
				if (rundays[1]<ftide1[i,1])
					break
				end
				fele1 = ftide1[i,2]+(ftide1[i+1,2]-ftide1[i,2])*(rundays[1]-ftide1[i,1])/(ftide1[i+1,1]-ftide1[i,1])
			end
		
			lt=length(ftide2)
			for i = 1:lt-1
				if (rundays[1]<ftide2[i,1])
					break
				end
				fele2 = ftide2[i,2]+(ftide2[i+1,2]-ftide2[i,2])*(rundays[1]-ftide2[i,1])/(ftide2[i+1,1]-ftide2[i,1])
			end
		
			lt=length(ftide3)
			for i = 1:lt-1
				if (rundays[1]<ftide3[i,1])
					break
				end
				fele3 = ftide3[i,2]+(ftide3[i+1,2]-ftide3[i,2])*(rundays[1]-ftide3[i,1])/(ftide3[i+1,1]-ftide3[i,1])
			end
		
			lt=length(ftide4)
			for i = 1:lt-1
				if (rundays[1]<ftide4[i,1])
					break
				end
				fele4 = ftide4[i,2]+(ftide4[i+1,2]-ftide4[i,2])*(rundays[1]-ftide4[i,1])/(ftide4[i+1,1]-ftide4[i,1])
			end
		
			# tidal forcing at the ocean boundaries
			dele = (fele2-fele1)/(79-29)
			fele = fele1
			for i=29:79
				q[i,2,1] = fele *rho0*g
				q[i,1,1] = q[i,2,1]
				fele = fele + dele
			end

			dele = (fele4-fele3)/(14-5)
			fele = fele3
			for j=5:14
				q[86,j,1] = fele *rho0*g
				q[87,j,1] = q[86,j,1]
				fele = fele + dele
			end
			
		# solving Poisson equation						
			Poisson(np,drdx,drdy,drdz,nx,ny,nz,dx,dy,dz,nt,dt,dtin,g,rho0,omega,epsilon,f,
			u,v,w,un,vn,wn,uintegral,vintegral,diffu,diffv,diffw,p,q,dq,ele,advx,advy,advz,
			du,dv,dw,ustar,vstar,wstar,qstar,depth,dry,land,atop,abot,ae,aw,an,as,atotal,ib)
		
		# write data
		@sync if ( mod(runtime[1], outinter) < dtin[1] )
				#WriteRho(nx,ny,nz,rho)
				@async WriteU(nx,ny,nz,u,umm)
				@async WriteV(nx,ny,nz,v,vmm)
				@async WriteW(nx,ny,nz,w)
				@async WriteQ(nx,ny,nz,q)
				#WriteT(nx,ny,nz,t)
				#WriteS(nx,ny,nz,s)
				#WriteA(nx,ny,nz,age)
			# screen monitoring
			@async println("Data output at time = ", runtime[1]/86400.0, "   ", "dt = ", dt, "dtin = ", dtin[1])
			@async println("dry, eta, depth[42,25] =>  ", dry[42,25,2], "   ", ele[42,25], "   ", q[42,25,1]/(g*rho0),"   ", depth[42,25])
			@async println("u[42,25,2],v[42,25,2], w[42,25,2] =>  ", u[42,25,2], "   ", v[42,25,2], "   ", w[42,25,2])
			#println("advx[6,9,2], advy[6,9,2] =>  ", advx[6,9,2], "   ", advy[6,9,2])
			#println("diffu[6,9,2], diffv[6,9,2] =>  ", diffu[6,9,2], "   ", diffv[6,9,2])
			#println("t[6,9,2], s[6,9,2], rho[6,9,2] =>  ", t[6,9,2], "   ", s[6,9,2], "   ", rho[6,9,2])
		end
		
		if ( rundays[1] >= periods )
			break
		end
		
	end
############### simulation loop end
	
	close(csv1)
	close(csv2)
	close(csv3)
	close(csv4)
	close(csv5)
	close(csv6)
	close(csv7)
	close(csv8)
		