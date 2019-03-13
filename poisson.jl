# calculate right-hand side of Poisson equation

	function Poisson(nx,ny,nz,dx,dy,dz,nt,dt,dtin,g,rho0,omega,epsilon,f,u,v,w,un,vn,wn,uintegral,vintegral,diffu,diffv,diffw,
							p,q,dq,ele,advx,advy,advz,du,dv,dw,ustar,vstar,wstar,qstar,depth,dry,land,atop,abot,ae,aw,an,as,atotal,ib)
		
		drdx = 0.0
		drdy = 0.0
		drdz = 0.0
		
		pressx = 0.0
		pressy = 0.0
		pressz = 0.0
		alpha = 0.0
		
		ex1 = 0.0
		ex2 = 0.0
		ex3 = 0.0
		
		# prediction of ustar, vstar, wstar
		@inbounds for k = 2:nz+1, j = 2:ny+1, i = 2:nx+1
					
					if ( dry[i,j,k] == 0 )
						drdx = 1.0/(rho0*dx[i,j])
						drdy = 1.0/(rho0*dy[i,j])
						drdz = 1.0/(rho0*dz[k])
						
						um = 0.25*( u[i,j,k] + u[i-1,j,k] + u[i,j+1,k] + u[i-1,j+1,k] )
						vm = 0.25*( v[i,j,k] + v[i+1,j,k] + v[i,j-1,k] + v[i+1,j-1,k] )
						
						pressx = -drdx*( q[i+1,j,k] - q[i,j,k] ) -drdx*( p[i+1,j,k] - p[i,j,k] )
						alpha = f[j]*dtin
						if ( dry[i+1,j,k]==0 ) du[i,j,k] = cos(alpha)*u[i,j,k] + sin(alpha)*vm + dtin*pressx + advx[i,j,k] + diffu[i,j,k] - u[i,j,k] end
						
						pressy = -drdy*( q[i,j+1,k] - q[i,j,k] ) -drdy*( p[i,j+1,k] - p[i,j,k] )
						alpha = 0.5*( f[j] + f[j+1] )*dtin
						if ( dry[i,j+1,k]==0 ) dv[i,j,k] = cos(alpha)*v[i,j,k] - sin(alpha)*um + dtin*pressy + advy[i,j,k] + diffv[i,j,k] - v[i,j,k] end
						
						pressz = -drdz*( q[i,j,k-1] - q[i,j,k] )
						if ( dry[i,j,k-1]==0 ) dw[i,j,k] = dtin*pressz + advz[i,j,k] + diffw[i,j,k] end
					end
				
					
				ustar[i,j,k] = u[i,j,k]
				vstar[i,j,k] = v[i,j,k]
				wstar[i,j,k] = w[i,j,k]
		end
		
######### start external mode(testing, not used)
		for nrep = 1:nt
			@inbounds for k = 2:nz+1, j = 2:ny+1, i = 2:nx+1
				ustar[i,j,k] = ustar[i,j,k] + du[i,j,k]/real(nt)
				vstar[i,j,k] = vstar[i,j,k] + dv[i,j,k]/real(nt)
				wstar[i,j,k] = wstar[i,j,k] + dw[i,j,k]/real(nt)
			end
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
			@inbounds for k = 2:nz+1, j = 2:ny+1, i = 2:nx+1
					ex1 = (ustar[i,j,k]-ustar[i-1,j,k])*dz[k]
					ex2 = (vstar[i,j,k]-vstar[i,j-1,k])*dx[i,j]*dz[k]/dy[i,j]
					ex3 = (wstar[i,j,k]-wstar[i,j,k+1])*dx[i,j]
					qstar[i,j,k] = -1.0*rho0/dt*( ex1 + ex2 + ex3 )
			end
			
			SorLoop(nx,ny,nz,dx,dy,dz,dt,rho0,dry,g,q,dq,atop,abot,ae,aw,as,an,atotal,
							qstar,ustar,vstar,wstar,u,v,w,un,vn,wn,uintegral,vintegral,omega,epsilon)
							
		end 
######## end external mode

			Update(nx,ny,nz,q,dq,u,v,w,un,vn,wn)
						# wet and dry update
			@. depth[1:nx+2,1:ny+2] = ele[1:nx+2,1:ny+2] + q[1:nx+2,1:ny+2,1]/(g*rho0)
			#WetDry(nx,ny,nz,dx,dy,dz,u,v,w,q,dq,depth,dry,land,atop,abot,ae,aw,an,as,atotal,ib)
			
    end
    # end function Poisson
    
    # sor(successive over relaxation) main function
	function SorLoop(nx,ny,nz,dx,dy,dz,dt,rho0,dry,g,q,dq,atop,abot,ae,aw,as,an,atotal,
							qstar,ustar,vstar,wstar,u,v,w,un,vn,wn,uintegral,vintegral,omega,epsilon)
			
		nmax = 8000	# maximum number of iterations
		# ntr = Threads.nthreads() # not used
		perror = zeros(Float64, nmax)
		
		@inbounds for n = 1:nmax
			
			#Threads.@threads for nn = 1:ntr
			
			perror[n] = 0.0
			
			Step1(n,nx,ny,nz,dry,dq,atop,abot,ae,aw,as,an,atotal,qstar,perror,omega)
			Step2(nx,ny,nz,dx,dy,dz,dt,dry,dq,un,vn,wn,ustar,vstar,wstar,rho0)
			Step3(nx,ny,nz,dz,uintegral,vintegral,un,vn)
			Step4(nx,ny,dx,dy,dq,dt,rho0,g,uintegral,vintegral)
							
			if ( perror[n] <= epsilon )
				#println("Number of Iterations =>", n, " Perror = ", perror[n])
				#nper = n
				break
			end
			
			
		end
		
		@label finish
		
		@inbounds for k = 2:nz+1
			for j = 2:ny+1
				for i = 2:nx+1
				ustar[i,j,k] = un[i,j,k]
				vstar[i,j,k] = vn[i,j,k]
				wstar[i,j,k] = wn[i,j,k]
				end
			end
		end
		
    end
    
    # Iteration
	function Step1(n,nx,ny,nz,dry,dq,atop,abot,ae,aw,as,an,atotal,qstar,perror,omega)
			
			q1 = 0.0
			q2 = 0.0
			term1 = 0.0
			term2 = 0.0
			term3 = 0.0
			term4 = 0.0
			
			# STEP 1 : predict delta q error
			@inbounds for k = 2:nz+1
				for j = 2:ny+1
					for i = 2:nx+1
						if ( dry[i,j,k] == 0 )
							q1 = dq[i,j,k]
							term1 = (atop[i,j,k] * dq[i,j,k-1]) + (abot[i,j,k] * dq[i,j,k+1])
							term2 = (aw[i,j,k] * dq[i-1,j,k]) + (ae[i,j,k] * dq[i+1,j,k])
							term3 = (as[i,j,k] * dq[i,j-1,k]) + (an[i,j,k] * dq[i,j+1,k])
							term4 = qstar[i,j,k] + term1 + term2 + term3
							q2 = (1.0-omega)*q1 + (omega*term4)/atotal[i,j,k]
							dq[i,j,k] = q2
							perror[n] = max( abs(q2-q1), perror[n])
						end
					end
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
	
	function Step2(nx,ny,nz,dx,dy,dz,dt,dry,dq,un,vn,wn,ustar,vstar,wstar,rho0)
			
			drdx = 0.0
			drdy = 0.0
			drdz = 0.0
			pressx = 0.0
			pressy = 0.0
			pressz = 0.0
		
			# STEP 2 : predict next timestep velocities
			@inbounds for k = 2:nz+1, j = 2:ny+1, i = 2:nx+1
			
						if ( dry[i,j,k] == 0 )
							drdx = 1.0/(rho0*dx[i,j])
							drdy = 1.0/(rho0*dy[i,j])
							drdz = 1.0/(rho0*dz[k])
							pressx = -drdx*( dq[i+1,j,k] - dq[i,j,k] )
							if ( dry[i+1,j,k] == 0 ) un[i,j,k] = ustar[i,j,k] + dt*pressx  end
							pressy = -drdy*( dq[i,j+1,k] - dq[i,j,k] )
							if ( dry[i,j+1,k] == 0 ) vn[i,j,k] = vstar[i,j,k] + dt*pressy  end
							pressz = -drdz*( dq[i,j,k-1] - dq[i,j,k] )
							if ( dry[i,j,k-1] == 0 ) wn[i,j,k] = wstar[i,j,k] + dt*pressz  end
						end
			end
	end
	
	function Step3(nx,ny,nz,dz,uintegral,vintegral,un,vn)
			
			# STEP 3 : predict vertical-intergrated flow(m^2/s)
			@inbounds for j = 2:ny+1, i = 2:nx+1
					uintegral[i,j] = 0.0
					vintegral[i,j] = 0.0
					
					@simd for k = 2:nz+1
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
			end
	end
	
	function Update(nx,ny,nz,q,dq,u,v,w,un,vn,wn)
	
		# update
		@inbounds for k = 2:nz+1, j = 2:ny+1, i = 2:nx+1
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
	
	function WetDry(nx,ny,nz,dx,dy,dz,u,v,w,q,dq,depth,dry,land,atop,abot,ae,aw,an,as,atotal,ib)
	@. dry = 0.0
	# wet and dry cell setup. if dry=0 is wet. else dry~=0 is dry.
	 	@inbounds for j = 1:ny+2, i = 1:nx+2
 				if (depth[i,j] <= 0.0) 
 					land[i,j] = 1
					nb = 0
 				else
 					zsum = 0.0
 					na = trunc(Int64,depth[i,j])
 					for k = 2:nz+1
 						nb = k-1
	 					zsum = zsum+dz[k]
	 					if ( na <= zsum )
		 					break
	 					end
	 				end
 					nb = min(nb,nz)
 				end
 				
				for k = nb+2:nz+2	# this point is very important nb+2, nz+2
					dry[i,j,k] = 1
					u[i,j,k] = 0.0
					v[i,j,k] = 0.0
					w[i,j,k] = 0.0
					q[i,j,k] = q[i,j,k] - dq[i,j,k]
				end
 				# dry(i,j,1) = 1	# rigid-lid approximation(not used)
	 	end
	
		# reset the poisson equation parameters					
		@inbounds for i = 2:nx+1, j = 2:ny+1, k = 2:nz+1
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
	
		# reset the ib(bottom index)
		@inbounds for k = 1:nz+1, j = 1:ny+2, i = 1:nx+2
					if (dry[i,j,k] == 0) && (dry[i,j,k+1] != 0) ib[i,j] = k end
		end

end