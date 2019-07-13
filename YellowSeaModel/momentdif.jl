# momentum diffusion.

	function MomentumDiffusion(nx,ny,nz,dx,dy,dz,dtin,ib,u,v,w,ah,az,dry,land,slip,taubu,diffu,taubv,diffv,diffw)
		
		dife = 0.0
		difw = 0.0
		difn = 0.0
		difs = 0.0
		dift = 0.0
		difb = 0.0
		difhor = 0.0
		difver = 0.0
		
		ahe = 0.0
		ahw = 0.0
		ahn = 0.0
		ahs = 0.0
		aztop = 0.0
		azbot = 0.0
			
		@fastmath @inbounds for k = 2:nz+1, j = 2:ny+1, i = 2:nx+1
					# u-diffusion
					ahe = ah[i+1,j,k]
					dife = ahe*( u[i+1,j,k] - u[i,j,k] )/dx[i,j]
					if ( dry[i+1,j,k] != 0 ) dife = 0.0 end
					ahw = ah[i,j,k]
					difw = ahw*( u[i,j,k] - u[i-1,j,k] )/dx[i,j]
					if ( dry[i-1,j,k] != 0 ) difw = 0.0 end
					
					ahn = 0.25*( ah[i,j,k] + ah[i+1,j,k] + ah[i,j+1,k] + ah[i+1,j+1,k] )
					difn = ahn*( u[i,j+1,k] - u[i,j,k] )/dy[i,j]
					if ( dry[i,j+1,k] !=0 ) difn = 0.0 end
					if ( land[i,j+1] !=0 ) difn = ahn*( u[i,j+1,k] - slip*u[i,j,k] )/dy[i,j] end
					
					ahs = 0.25*( ah[i,j,k] + ah[i+1,j,k] + ah[i,j-1,k] + ah[i+1,j-1,k] )
					difs = ahs*( u[i,j,k] - u[i,j-1,k] )/dy[i,j]
					if ( dry[i,j-1,k] !=0 ) difs = 0.0 end
					if ( land[i,j-1] !=0 ) difs = ahs*( slip*u[i,j,k] - u[i,j-1,k] )/dy[i,j] end
					
					difhor = ( dife-difw )/dx[i,j] + ( difn-difs )/dy[i,j]
					
					aztop = 0.25*( az[i,j,k] + az[i,j,k-1] + az[i+1,j,k] + az[i+1,j,k-1] )
					dift = aztop*( u[i,j,k-1] - u[i,j,k] )/dz[k]
					azbot = 0.25*( az[i,j,k] + az[i,j,k+1] + az[i+1,j,k] + az[i+1,j,k+1] )
					difb = azbot*( u[i,j,k] - u[i,j,k+1] )/dz[k]
					if ( k == ib[i,j] ) difb = taubu[i,j] end
					difver = ( dift - difb )/dz[k]
					
					diffu[i,j,k] =  dtin[1]*( difhor + difver )
					
					# v-diffusion
					ahe = 0.25*( ah[i,j,k] + ah[i,j+1,k] + ah[i+1,j,k] + ah[i+1,j+1,k] )
					dife = ahe*( v[i+1,j,k] - v[i,j,k] )/dx[i,j]
					if ( dry[i+1,j,k] !=0 ) dife = 0.0 end
					if ( land[i+1,j] !=0 ) dife = ahe*( v[i+1,j,k] - slip*v[i,j,k] )/dx[i,j] end
					
					ahw = 0.25*( ah[i,j,k] + ah[i,j+1,k] + ah[i-1,j,k] + ah[i-1,j+1,k] )
					difw = ahw*( v[i,j,k] - v[i-1,j,k] )/dx[i,j]
					if ( dry[i-1,j,k] !=0 ) difw = 0.0 end
					if ( land[i-1,j] !=0 ) difw = ahw*( slip*v[i,j,k] - v[i-1,j,k] )/dx[i,j] end
					
					ahn = ah[i,j+1,k]
					difn = ahn*( v[i,j+1,k] - v[i,j,k] )/dy[i,j]
					if ( dry[i,j+1,k] != 0 ) difn = 0.0 end
					ahs = ah[i,j,k]
					difs = ahs*( v[i,j,k] - v[i,j-1,k] )/dy[i,j]	#
					if ( dry[i,j-1,k] != 0 ) difs = 0.0 end
					
					difhor = ( dife-difw )/dx[i,j] + ( difn-difs )/dy[i,j]
					
					aztop = 0.25*( az[i,j,k] + az[i,j+1,k] + az[i,j,k-1] + az[i,j+1,k-1] )
					dift = aztop*( v[i,j,k-1] - v[i,j,k] )/dz[k]
					azbot = 0.25*( az[i,j,k] + az[i,j+1,k] + az[i,j,k+1] + az[i,j+1,k+1] )
					difb = azbot*( v[i,j,k] - v[i,j,k+1] )/dz[k]
					if ( k == ib[i,j] ) difb = taubv[i,j] end
					difver = ( dift - difb )/dz[k]
					
					diffv[i,j,k] =  dtin[1]*( difhor + difver )
					
					# w-diffusion
					ahe = 0.25*( ah[i,j,k] + ah[i,j,k-1] + ah[i+1,j,k] + ah[i+1,j,k-1] )
					dife = ahe*( w[i+1,j,k] - w[i,j,k] )/dx[i,j]
					if ( dry[i+1,j,k] !=0 ) dife = 0.0 end

					ahw = 0.25*( ah[i,j,k] + ah[i,j,k-1] + ah[i-1,j,k] + ah[i-1,j,k-1] )
					difw = ahw*( w[i,j,k] - w[i-1,j,k] )/dx[i,j]
					if ( dry[i-1,j,k] !=0 ) difw = 0.0 end
					
					ahn = 0.25*( ah[i,j,k] + ah[i,j,k-1] + ah[i,j+1,k] + ah[i,j+1,k-1] )
					difn = ahn*( w[i,j+1,k] - w[i,j,k] )/dy[i,j]
					if ( dry[i,j+1,k] != 0 ) difn = 0.0 end
					
					ahs = 0.25*( ah[i,j,k] + ah[i,j,k-1] + ah[i,j-1,k] + ah[i,j-1,k-1] )
					difs = ahs*( w[i,j,k] - w[i,j-1,k] )/dy[i,j]
					if ( dry[i,j-1,k] != 0 ) difs = 0.0 end
					
					difhor = ( dife-difw )/dx[i,j] + ( difn-difs )/dy[i,j]
					
					aztop = az[i,j,k-1]
					dift = aztop*( w[i,j,k-1] - w[i,j,k] )/dz[k]
					if ( dry[i,j,k-1] != 0 ) dift = 0.0 end
					azbot = az[i,j,k]
					difb = azbot*( w[i,j,k] - w[i,j,k+1] )/dz[k]
					if ( dry[i,j,k+1] != 0 ) difb = 0.0 end
					difver = ( dift - difb )/dz[k]
					
					diffw[i,j,k] =  dtin[1]*( difhor + difver )
		end
		
    end
    
  
