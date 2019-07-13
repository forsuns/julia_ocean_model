# calculate right-hand side of Poisson equation

	function Poisson(np,drdx,drdy,drdz,nx,ny,nz,dx,dy,dz,nt,dt,dtin,g,rho0,omega,epsilon,f,u,v,w,un,vn,wn,uintegral,vintegral,diffu,diffv,diffw,
						p,q,dq,ele,advx,advy,advz,du,dv,dw,ustar,vstar,wstar,qstar,depth,dry,land,atop,abot,ae,aw,an,as,atotal,ib)
		
		pressx = 0.0
		pressy = 0.0
		pressz = 0.0
		alpha = 0.0
		
		ex1 = 0.0
		ex2 = 0.0
		ex3 = 0.0
		
		um = 0.0
		vm = 0.0
		
		# prediction of ustar, vstar, wstar
		@fastmath @inbounds for j = 2:ny+1, i = 2:nx+1, k = 2:ib[i,j]
					if ( dry[i,j,k] == 0 )
						um = 0.25*( u[i,j,k] + u[i-1,j,k] + u[i,j+1,k] + u[i-1,j+1,k] )
						vm = 0.25*( v[i,j,k] + v[i+1,j,k] + v[i,j-1,k] + v[i+1,j-1,k] )
						
						pressx = -drdx[i,j]*( q[i+1,j,k] - q[i,j,k] ) -drdx[i,j]*( p[i+1,j,k] - p[i,j,k] )
						alpha = f[j]* dtin[1]
						if ( dry[i+1,j,k]==0 ) du[i,j,k] = cos(alpha)*u[i,j,k] + sin(alpha)*vm +  dtin[1]*pressx + advx[i,j,k] + diffu[i,j,k] - u[i,j,k] end
						
						pressy = -drdy[i,j]*( q[i,j+1,k] - q[i,j,k] ) -drdy[i,j]*( p[i,j+1,k] - p[i,j,k] )
						alpha = 0.5*( f[j] + f[j+1] )* dtin[1]
						if ( dry[i,j+1,k]==0 ) dv[i,j,k] = cos(alpha)*v[i,j,k] - sin(alpha)*um +  dtin[1]*pressy + advy[i,j,k] + diffv[i,j,k] - v[i,j,k] end
						
						pressz = -drdz[k]*( q[i,j,k-1] - q[i,j,k] )
						if ( dry[i,j,k-1]==0 ) dw[i,j,k] =  dtin[1]*pressz + advz[i,j,k] + diffw[i,j,k] end
					end	
				ustar[i,j,k] = u[i,j,k]
				vstar[i,j,k] = v[i,j,k]
				wstar[i,j,k] = w[i,j,k]
		end
		
		np[1] = 0
		
			External(nt,np,nx,ny,nz,dx,dy,dz,dt,rho0,dry,g,q,dq,atop,abot,ae,aw,as,an,atotal,drdx,drdy,drdz,
				qstar,ustar,vstar,wstar,u,v,w,un,vn,wn,uintegral,vintegral,omega,epsilon,ib)

			Update(nx,ny,nz,q,dq,u,v,w,un,vn,wn)
			
    end
    # end function Poisson
    
    ######### start external mode function(testing, not used)
    function External(nt,np,nx,ny,nz,dx,dy,dz,dt,rho0,dry,g,q,dq,atop,abot,ae,aw,as,an,atotal,drdx,drdy,drdz,
				qstar,ustar,vstar,wstar,u,v,w,un,vn,wn,uintegral,vintegral,omega,epsilon,ib)
    
		for nrep = 1:nt[1]
			@. ustar = ustar+du/real(nt[1])
			@. vstar = vstar+dv/real(nt[1])
			@. wstar = wstar+dw/real(nt[1])
			# boundary
			@inbounds for k = 2:nz+1, j = 1:ny+2
				ustar[1,j,k] = 0.0
				ustar[nx+1,j,k] = 0.0
				vstar[1,j,k] = vstar[2,j,k]
				vstar[nx+2,j,k] = vstar[nx+1,j,k]
				wstar[1,j,k] = wstar[2,j,k]
				wstar[nx+2,j,k] = wstar[nx+1,j,k]
			end
		
			@inbounds for k = 2:nz+1, i = 1:nx+2
				ustar[i,1,k] = ustar[i,2,k]
				ustar[i,ny+2,k] = ustar[i,ny+1,k]
				vstar[i,1,k] = vstar[i,2,k]
				vstar[i,ny+1,k] = vstar[i,ny,k]
				wstar[i,1,k] = wstar[i,2,k]
				wstar[i,ny+2,k] = wstar[i,ny+1,k]
			end
		
			# prediction of qstar
			@fastmath @inbounds for j = 2:ny+1, i = 2:nx+1, k = 2:ib[i,j]
							ex1 = (ustar[i,j,k]-ustar[i-1,j,k])*dz[k]
							ex2 = (vstar[i,j,k]-vstar[i,j-1,k])*dx[i,j]*dz[k]/dy[i,j]
							ex3 = (wstar[i,j,k]-wstar[i,j,k+1])*dx[i,j]
							qstar[i,j,k] = -1.0*rho0/dt*( ex1 + ex2 + ex3 )
			end
			
			SorLoop(np,nx,ny,nz,dx,dy,dz,dt,rho0,dry,g,q,dq,atop,abot,ae,aw,as,an,atotal,drdx,drdy,drdz,
							qstar,ustar,vstar,wstar,u,v,w,un,vn,wn,uintegral,vintegral,omega,epsilon,ib)
		end 
	end
	######## end external mode function
	
    # sor(successive over relaxation) main function
	function SorLoop(np,nx,ny,nz,dx,dy,dz,dt,rho0,dry,g,q,dq,atop,abot,ae,aw,as,an,atotal,drdx,drdy,drdz,
							qstar,ustar,vstar,wstar,u,v,w,un,vn,wn,uintegral,vintegral,omega,epsilon,ib)
							
		nmax =8000	# maximum number of iterations
		# ntr = Threads.nthreads() # not used
		perror = zeros(Float64, 1)
		
		@inbounds for n = 1:nmax
			#Threads.@threads for nn = 1:ntr
			
			perror[1] = 0.0
			
			Step1(n,nx,ny,nz,dry,land,dq,atop,abot,ae,aw,as,an,atotal,qstar,perror,omega,ib)
			Step2(nx,ny,nz,dx,dy,dz,dt,dry,land,dq,un,vn,wn,ustar,vstar,wstar,rho0,drdx,drdy,drdz,ib)
			Step3(nx,ny,nz,dz,uintegral,vintegral,un,vn,ib)
			Step4(nx,ny,dx,dy,dq,dt,rho0,g,uintegral,vintegral)
			
			if ( perror[1] <= epsilon )
				#println("Number of Iterations =>", n, " Perror = ", perror[1])
				np[1] = max(np[1],n)
				#@goto finish
				break
			end
			
			if ( n == nmax)	println("Arrived the maximum iteration number ",n," error ",perror[1]) end
			
		end
		
		#@label finish
		
		@fastmath @inbounds for j = 2:ny+1, i = 2:nx+1, k = 2:ib[i,j]
						ustar[i,j,k] = un[i,j,k]
						vstar[i,j,k] = vn[i,j,k]
						wstar[i,j,k] = wn[i,j,k]
		end
		
    end
	
    # Iteration
	function Step1(n,nx,ny,nz,dry,land,dq,atop,abot,ae,aw,as,an,atotal,qstar,perror,omega,ib)
			
			q1 = 0.0
			q2 = 0.0
			term1 = 0.0
			term2 = 0.0
			term3 = 0.0
			term4 = 0.0
			
			# STEP 1 : predict delta q error
			@fastmath @inbounds for j = 2:ny+1, i = 2:nx+1, k = 2:ib[i,j]
								if ( dry[i,j,k] == 0 )
								q1 = dq[i,j,k]
								term1 = (atop[i,j,k] * dq[i,j,k-1]) + (abot[i,j,k] * dq[i,j,k+1])
								term2 = (aw[i,j,k] * dq[i-1,j,k]) + (ae[i,j,k] * dq[i+1,j,k])
								term3 = (as[i,j,k] * dq[i,j-1,k]) + (an[i,j,k] * dq[i,j+1,k])
								term4 = qstar[i,j,k] + term1 + term2 + term3
								q2 = (1.0-omega)*q1 + (omega*term4)/atotal[i,j,k]
								dq[i,j,k] = q2
								perror[1] = max( abs(q2-q1), perror[1])
								end
			end
			
			@inbounds for k = 2:nz+1, j = 1:ny+2
						dq[1,j,k] = dq[2,j,k]
						dq[nx+2,j,k] = dq[nx+1,j,k]
			end
			@inbounds for k = 2:nz+1, i = 1:nx+2
						dq[i,1,k] = dq[i,2,k]
						dq[i,ny+2,k] = dq[i,ny+1,k]
			end
			
	end
	
	function Step2(nx,ny,nz,dx,dy,dz,dt,dry,land,dq,un,vn,wn,ustar,vstar,wstar,rho0,drdx,drdy,drdz,ib)
			
			pressx = 0.0
			pressy = 0.0
			pressz = 0.0
		
			# STEP 2 : predict next timestep velocities
			@fastmath @inbounds for j = 2:ny+1, i = 2:nx+1, k = 2:ib[i,j]
							un[i,j,k] = 0.0
							vn[i,j,k] = 0.0
							wn[i,j,k] = 0.0
							if ( dry[i,j,k] == 0 )
								pressx = -drdx[i,j]*( dq[i+1,j,k] - dq[i,j,k] )
								if ( dry[i+1,j,k] == 0 ) un[i,j,k] = ustar[i,j,k] + dt*pressx  end
								pressy = -drdy[i,j]*( dq[i,j+1,k] - dq[i,j,k] )
								if ( dry[i,j+1,k] == 0 ) vn[i,j,k] = vstar[i,j,k] + dt*pressy  end
								pressz = -drdz[k]*( dq[i,j,k-1] - dq[i,j,k] )
								if ( dry[i,j,k-1] == 0 ) wn[i,j,k] = wstar[i,j,k] + dt*pressz  end
							end
			end
	end
	
	function Step3(nx,ny,nz,dz,uintegral,vintegral,un,vn,ib)
			
			# STEP 3 : predict vertical-intergrated flow(m^2/s)
			@inbounds for j = 2:ny+1, i = 2:nx+1
					uintegral[i,j] = 0.0
					vintegral[i,j] = 0.0
					
					@simd for k = 2:ib[i,j]
						uintegral[i,j] = uintegral[i,j] + dz[k]*un[i,j,k]
						vintegral[i,j] = vintegral[i,j] + dz[k]*vn[i,j,k]
					end
			end
			
			# lateral boundary
			@inbounds for j = 2:ny+1
				uintegral[1,j] = 0.0
				uintegral[nx+1,j] = 0.0
			end
			
			@inbounds for i = 2:nx+1
				vintegral[i,ny+1] = vintegral[i,ny]
				vintegral[i,1] = vintegral[i,2]
			end
			
	end
	
	
	function Step4(nx,ny,dx,dy,dq,dt,rho0,g,uintegral,vintegral)
			
			# STEP 4 : predict dynamic pressure field at surface
			@inbounds for j = 2:ny+1, i = 2:nx+1
					dq[i,j,1] = -dt*rho0*g*( (uintegral[i,j] - uintegral[i-1,j])/dx[i,j] + ( vintegral[i,j] - vintegral[i,j-1])/dy[i,j] )
					if ( dry[i,j,2] != 0.0 ) dq[i,j,1] = 0.0 end
			end
	end
	
	function Update(nx,ny,nz,q,dq,u,v,w,un,vn,wn)
	
		# update
		@fastmath @inbounds for j = 2:ny+1, i = 2:nx+1, k = 2:ib[i,j]
						q[i,j,k] = q[i,j,k] + dq[i,j,k]
						u[i,j,k] = un[i,j,k]
						v[i,j,k] = vn[i,j,k]
						w[i,j,k] = wn[i,j,k]
		end
		
		@inbounds for j = 2:ny+1, i = 2:nx+1
				q[i,j,1] = q[i,j,1] + dq[i,j,1]
		end
		
		# lateral boundary
		@inbounds for j = 2:ny+1
			q[1,j,1] = q[2,j,1]
			q[nx+2,j,1] = q[nx+1,j,1]
			@simd for k = 2:nz+1
				u[1,j,k] = 0.0
				u[nx+1,j,k] = 0.0
				v[1,j,k] = v[2,j,k]
				v[nx+2,j,k] = v[nx+1,j,k]
				w[1,j,k] = w[2,j,k]
				w[nx+2,j,k] = w[nx+1,j,k]
				q[1,j,k] = q[2,j,k]
				q[nx+2,j,k] = q[nx+1,j,k]
			end
		end
		
		@inbounds for i = 2:nx+1
			q[i,1,1] = q[i,2,1]
			q[i,ny+2,1] = q[i,ny+1,1]
			@simd for k = 2:nz+1
				u[i,1,k] = u[i,2,k]
				u[i,ny+2,k] = u[i,ny+1,k]
				v[i,1,k] = v[i,2,k]
				v[i,ny+1,k] = v[i,ny,k]
				w[i,1,k] = w[i,2,k]
				w[i,ny+2,k] = w[i,ny+1,k]
				q[i,1,k] = q[i,2,k]
				q[i,ny+2,k] = q[i,ny+1,k]
			end
		end
		
	end
	
	function WetDry(nx,ny,nz,u,v,w,q,depth)
	# wet and dry cell setup. if dry=0 is wet. else dry~=0 is dry.
	 	@inbounds for j = 2:ny+1, i = 2:nx+1
	 		if ( depth[i,j] <= 0.0 )
	 			for k = 2:nz+1
	 				u[i,j,k] = 0.0
	 				v[i,j,k] = 0.0
	 				w[i,j,k] = 0.0
	 				q[i,j,k] = 0.0
	 			end
	 		end
	 	end

	end