# eddy vertical

	function EddyVertical(nx,ny,nz,dz,u,v,rho,kz,az,rho0,g)
	
		vismin = 0.0001
		vismax = 0.05
		windmix = vismax
		c1 = 0.15
		c2 = 0.0
		utop = 0.0
		ubot = 0.0
		dudz = 0.0
		eddyx = 0.0
		vtop = 0.0
		vbot = 0.0
		dvdz = 0.0
		eddyy = 0.0
		stab = 0.0
		wurz = 0.0
              
        @inbounds for i = 2:nx+1, j = 2:ny+1, k = 2:nz+1
        		
	        		c2 = c1*c1*dz[k]*dz[k]
        			utop = 0.5 * ( u[i,j,k-1] + u[i-1,j,k-1] )
        			ubot = 0.5 * ( u[i,j,k+1] + u[i-1,j,k+1] )
        			dudz = (utop-ubot)/(2.0*dz[k])
        			eddyx = dudz*dudz
        			
        			vtop = 0.5 * ( v[i,j,k-1] + v[i,j-1,k-1] )
        			vbot = 0.5 * ( v[i,j,k+1] + v[i,j-1,k+1] )
        			dvdz = (vtop-vbot)/(2.0*dz[k])
        			eddyy = dvdz*dvdz
        			
        			stab = -g/rho0*(rho[i,j,k-1]-rho[i,j,k+1])/(2.0*dz[k])
        			wurz = max(0.0, eddyx+eddyy-stab)
        			wurz = max(c2*sqrt(wurz),vismin)
        			wurz = min(wurz,vismax)
        			if (stab<0) wurz = vismax end
        			if (k<3) wurz = windmix end
        			kz[i,j,k] = wurz
        end
        
        @inbounds for j = 2:ny+1, i = 2:nx+1
        		kz[i,j,1] = kz[i,j,2]
        		kz[i,j,nz+2] = kz[i,j,nz+1]
        end
        
        @inbounds for k = 1:nz+2, i = 2:nx+1
        		kz[i,1,k] = kz[i,2,k]
        		kz[i,ny+2,k] = kz[i,ny+1,k]
        end
        
        @inbounds for k = 1:nz+2, j = 2:ny+1
        		kz[1,j,k] = kz[2,j,k]
        		kz[nx+2,j,k] = kz[nx+1,j,k]
        end
        
		@. az = kz
        
    end