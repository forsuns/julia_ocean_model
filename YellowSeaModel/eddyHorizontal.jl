# Smagorinsky(1963) eddy horizontal

	function EddyHorizontal(nx,ny,nz,dx,dy,u,v,w,ah,kh)
	
		c1 = 0.1
		c2 = 0.0
		term1 = 0.0
		term2 = 0.0
		term3 = 0.0
		uNorth = 0.0
		uSouth = 0.0
		vEast = 0.0
		vWest = 0.0
              
        @fastmath @inbounds for i = 2:nx+1, j = 2:ny+1, k = 2:nz+1
        		
   				c2 = c1*dx[i,j]*dy[i,j]
        		term1 = (u[i,j,k] - u[i-1,j,k])/dx[i,j]
        		term1 = term1 * term1
        		
        		term2 = (v[i,j,k] - v[i,j-1,k])/dy[i,j]
        		term2 = term2 * term2
        		
        		uNorth = 0.25*(u[i,j,k]+u[i-1,j,k]+u[i,j+1,k]+u[i-1,j+1,k])
				uSouth = 0.25*(u[i,j,k]+u[i-1,j,k]+u[i,j-1,k]+u[i-1,j-1,k])
				vEast = 0.25*(v[i,j,k]+v[i,j-1,k]+v[i+1,j,k]+v[i+1,j-1,k])
				vWest = 0.25*(v[i,j,k]+v[i,j-1,k]+v[i-1,j,k]+v[i-1,j-1,k])
				
				term3 = (uNorth-uSouth)/dy[i,j] + (vEast-vWest)/dx[i,j]
				term3 = 0.5*term3*term3
				
				ah[i,j,k] = c2*sqrt(term1+term2+term3) + 5.0
        end
        
        @inbounds for j = 2:ny+1, i = 2:nx+1
        		ah[i,j,1] = ah[i,j,2]
        		ah[i,j,nz+2] = ah[i,j,nz+1]
        end
        
        @inbounds for k = 1:nz+2, i = 2:nx+1
        		ah[i,1,k] = ah[i,2,k]
        		ah[i,ny+2,k] = ah[i,ny+1,k]
        end
        
        @inbounds for k = 1:nz+2, j = 2:ny+1
        		ah[1,j,k] = ah[2,j,k]
        		ah[nx+2,j,k] = ah[nx+1,j,k]
        end
        
		@. kh = ah
        
    end