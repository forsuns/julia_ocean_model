# main program

include("readtopo.jl")
include("eddyVertical.jl")
include("eddyHorizontal.jl")
include("advection.jl")
include("density.jl")
include("hydrostatic.jl")
include("momentadv.jl")
include("bottomShear.jl")
include("surfaceWindShear.jl")
include("momentdif.jl")
include("poisson.jl")

# for parallel computing(not used)
#@everywhere using Distributed
#@everywhere using DistributedArrays
#@everywhere using DistributedArrays.SPMD

# initialize

	nx = 9
	ny = 15
	nz = 10
	dt = 2.0	# exterior time step
	dtin = dt	# interior time step(not used)
	nt = trunc(Int64,dtin/dt)
	g = 9.81			# gravity acceleration
	rho0 = 1028.0		# background rho
	warmup = 0.0		# warm up parameter
	drag = 0.100		# bottom friction rate
	slip = 1.0			# slip-nonslip condition indicator
	fbeta = 2.2e-11		# coriolis parameter
	omega = 1.4			# S.O.R parameter
	epsilon = 0.0001			# possion accuracy parameter

	periods = 1.0				# total running time (days)
	outinter = 1.0/6.0			# output interval (hours)
	runtime = 0.0
	nper = 1				# iteration default number
	
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
	um = zeros(Float64, nz+1)
	vm = zeros(Float64, nz+1)
	
	# Eddy Viscosity Variables
	kz = zeros(Float64, nx+2,ny+2,nz+2)
	az = zeros(Float64, nx+2,ny+2,nz+2)
	kh = zeros(Float64, nx+2,ny+2,nz+2)
	ah = zeros(Float64, nx+2,ny+2,nz+2)
	
	# Advection Variables
	rho = zeros(Float64, nx+2,ny+2,nz+2)
	rhon = zeros(Float64, nx+2,ny+2,nz+2)
	B = zeros(Float64, nx+2,ny+2,nz+2)
	BN = zeros(Float64, nx+2,ny+2,nz+2)
	#B = dzeros((nx+2,ny+2,nz+2),workers()[1:4],[1,1,4])
	#BN = dzeros((nx+2,ny+2,nz+2),workers()[1:4],[1,1,4])
	
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
	
		@. dx = 100.0
		@. dy = 100.0
		@. dz = 1.0
		@. ah = 1.0
		@. kh = 1.0

		# coriolis parameter setup
		@. f = 0.0000
		
		# read topography from file
		ReadTopo(ele,ib)
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
		# open iterations number n.csv
		fname = "./output/n.csv"
		csv5 = open(fname,"w+")
		# open um.csv
		fname = "./output/um.csv"
		csv6 = open(fname,"w+")
		# open vm.csv
		fname = "./output/vm.csv"
		csv7 = open(fname,"w+")
		# open q.csv
		fname = "./output/q.csv"
		csv8 = open(fname,"w+")
		
		write(csv1, ["i, j, rho_1 ~ rho_k \n"])
		write(csv2, ["i, j, u_1 ~ u_k \n"])
		write(csv3, ["i, j, v_1 ~ v_k \n"])
		write(csv4, ["i, j, w_1 ~ w_k \n"])
		write(csv8, ["i, j, q_1 ~ q_k \n"])
		
		@simd for j = 1:ny+2
			ele[nx+2,j] = 0.0
			ele[1,j] = 0.0
		end
		
		@simd for i = 1:nx+2
			ele[i,ny+2] = 0.0
			ele[i,1] = 0.0
		end
	
		# initializing wet and dry
		@. depth = ele
		WetDry(nx,ny,nz,dx,dy,dz,u,v,w,q,dq,depth,dry,land,atop,abot,ae,aw,an,as,atotal,ib)
		
	# wind-stress forcing
	taux .= 0.0
	tauy .= 0.0
			
	# runtime set
	ntotal = trunc(Int64,periods*24*60*60/dtin)
	runtime = 0.0
	nout = trunc(Int64,outinter*60*60/dtin)
		
	# output initial values	
	# write h.csv
		fname = "./output/h.csv"
		csvfile = open(fname,"w")
		write(csvfile, ["i, j, depth, land, ib\n"])
		@inbounds for j = 1:ny+2, i = 1:nx+2
				write(csvfile, [string.(i),",",string.(j),",",string.(depth[i,j]),",",string.(land[i,j]),",",string.(ib[i,j]),"\n"])
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
					write(csv7, [string.(i),",",string.(j),",",string.(q[i,j,2:nz+1],","),"\n"])
			end
	
############### simulation loop start
	@inbounds for n = 1:ntotal
	
		runtime = runtime + dtin
		warmup = min( runtime/(1.0*24*3600), 1.0 )

		  q[nx+1,ny+1,1] = q[nx+1,ny+1,1] + ( (10.0 *dt ) / (dx[nx+1,ny+1]*dy[nx+1,ny+1]) ) *rho0*g
		  #q[nx+1,ny+2,1] = q[nx+1,ny+1,1]
		  #q[nx+2,ny+1,1] = q[nx+1,ny+1,1]
		  q[2,2,1] = q[2,2,1] - ( (10.0 *dt ) / (dx[2,2]*dy[2,2]) ) *rho0*g
		  #q[2,1,1] = q[2,2,1]
		  #q[1,2,1] = q[2,2,1]
		# flow dynamic calculate
			EddyVertical(nx,ny,nz,dz,u,v,rho,kz,az,rho0,g) 
			EddyHorizontal(nx,ny,nz,dx,dy,u,v,w,ah,kh) 
			BottomShear(nx,ny,ib,drag,u,v,taubu,taubv) 
			Hydrostatic(nx,ny,nz,dz,g,p,rho) 
			SurfaceWindShear(nx,ny,dz,warmup,rho0,az,u,v,taux,tauy) 
			Density(nx,ny,nz,dry,dx,dy,dz,dt,dtin,rho,rhon,B,BN,u,v,w,ah,kz) 
			MomentumAdvection(nx,ny,nz,dx,dy,dz,dt,dtin,u,v,w,B,BN,advx,advy,advz) 
			MomentumDiffusion(nx,ny,nz,dx,dy,dz,dt,dtin,ib,u,v,w,ah,az,dry,land,slip,taubu,diffu,taubv,diffv,diffw) 
			
		# tidal forcing at the ocean boundaries
		#	q[2:nx+1,ny+1,1] = -( 0.1*sin( 2.0*pi*(runtime+600.0)/(12.0*3600)) ) *rho0*g
		#	q[2:nx+1,ny+2,1] = q[2:nx+1,ny+1,1]
		#	q[2:nx+1,2,1] = -( 0.1*sin( 2.0*pi*(runtime+0.0)/(12.0*3600)) ) *rho0*g
		#	q[2:nx+1,1,1] = q[2:nx+1,2,1]
		
		# solving Poisson equation						
			Poisson(nx,ny,nz,dx,dy,dz,nt,dt,dtin,g,rho0,omega,epsilon,f,u,v,w,un,vn,wn,uintegral,vintegral,diffu,diffv,diffw,p,q,dq,ele,advx,advy,advz,du,dv,dw,ustar,vstar,wstar,qstar,depth,dry,land,atop,abot,ae,aw,an,as,atotal,ib)
		
		# write data
		if ( mod(n, nout) == 0 )
			
			@inbounds for j = 2:ny+1, i = 2:nx+1
					write(csv1, [string.(i),",",string.(j),",",string.(rho[i,j,2:nz+1],","),"\n"])
			end
			
			@inbounds for j = 2:ny+1, i = 2:nx+1
				um[2:nz+1] = 0.5*(u[i,j,2:nz+1] + u[i-1,j,2:nz+1])
				vm[2:nz+1] = 0.5*(v[i,j,2:nz+1] + v[i,j-1,2:nz+1])
					write(csv2, [string.(i),",",string.(j),",",string.(um[2:nz+1],","),"\n"])
					write(csv3, [string.(i),",",string.(j),",",string.(vm[2:nz+1],","),"\n"])
					write(csv4, [string.(i),",",string.(j),",",string.(w[i,j,2:nz+1],","),"\n"])
			end
			
			for j = 2:ny+1, i = 2:nx+1
					write(csv8, [string.(i),",",string.(j),",",string.(q[i,j,2:nz+1],","),"\n"])
			end
			
			# screen monitoring
			println("Data output at time = ", runtime/86400.0 )
			println("ele, eta, depth[6,9] =>  ", ele[6,9], "   ", q[6,9,1]/(g*rho0),"   ", depth[6,9])
			println("u[6,9,2], v[6,9,2] =>  ", u[6,9,2], "   ", v[6,9,2])
		end
		
	end
############### simulation loop start
	
	close(csv1)
	close(csv2)
	close(csv3)
	close(csv4)
	close(csv5)
	close(csv6)
	close(csv8)
		