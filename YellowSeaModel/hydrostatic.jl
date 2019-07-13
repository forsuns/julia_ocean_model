# hydrostatic pressure

	function Hydrostatic(nx,ny,nz,dz,g,p,rho)
		
		@inbounds for j = 1:ny+2, i = 1:nx+2
				p[i,j,1] = 0.0
				for k = 2:nz+2
					p[i,j,k] = p[i,j,k-1] + 0.5*( rho[i,j,k-1]+rho[i,j,k] )*g*dz[k]
				end
		end
    end