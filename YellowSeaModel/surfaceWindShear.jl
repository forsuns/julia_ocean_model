# surface wind shear stress

	function SurfaceWindShear(nx,ny,dz,warmup,rho0,az,u,v,taux,tauy)
		
		aztop = 0.0

		@fastmath @inbounds for j = 2:ny+1, i = 2:nx+1
				aztop = 0.25*( az[i,j,2]+az[i+1,j,2]+az[i,j,1]+az[i+1,j,1] )
				u[i,j,1] = u[i,j,2] + warmup*dz[2]*taux[i,j] / (rho0*aztop)
				aztop = 0.25*( az[i,j,2]+az[i,j+1,2]+az[i,j,1]+az[i,j+1,1] )
				v[i,j,1] = v[i,j,2] + warmup*dz[2]*tauy[i,j] / (rho0*aztop)
		end
		
    end