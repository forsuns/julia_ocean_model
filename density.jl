# density advection

	function Density(nx,ny,nz,dry,dx,dy,dz,dt,dtin,rho,rhon,B,BN,u,v,w,ah,kz)
		
		CuPos = zeros(Float64, nx+2,ny+2,nz+2)
		CuNeg = zeros(Float64, nx+2,ny+2,nz+2)
		CvPos = zeros(Float64, nx+2,ny+2,nz+2)
		CvNeg = zeros(Float64, nx+2,ny+2,nz+2)
		CwPos = zeros(Float64, nx+2,ny+2,nz+2)
		CwNeg = zeros(Float64, nx+2,ny+2,nz+2)
		
		divu = 0.0
		divv = 0.0
		divw = 0.0
		div = 0.0
		
		khe = 0.0
		khw = 0.0
		khn = 0.0
		khs = 0.0
		dife = 0.0
		difw = 0.0
		difn = 0.0
		difs = 0.0
		dift = 0.0
		difb = 0.0
		difhor = 0.0
		difver = 0.0
		dif = 0.0
		
		aztop = 0.0
		azbot = 0.0
		
		@inbounds for i = 1:nx+2, j = 1:ny+2, k = 1:nz+2
					CuPos[i,j,k] = 0.5*( u[i,j,k] + abs(u[i,j,k]) ) * dtin/dx[i,j]
					CuNeg[i,j,k] = 0.5*( u[i,j,k] - abs(u[i,j,k]) ) * dtin/dx[i,j]
					CvPos[i,j,k] = 0.5*( v[i,j,k] + abs(v[i,j,k]) ) * dtin/dy[i,j]
					CvNeg[i,j,k] = 0.5*( v[i,j,k] - abs(v[i,j,k]) ) * dtin/dy[i,j]
					CwPos[i,j,k] = 0.5*( w[i,j,k] + abs(w[i,j,k]) ) * dtin/dz[k]
					CwNeg[i,j,k] = 0.5*( w[i,j,k] - abs(w[i,j,k]) ) * dtin/dz[k]
		end
		
		@. B = rho
		
		Advection(nx,ny,nz,B,BN,CuPos,CuNeg,CvPos,CvNeg,CwPos,CwNeg)
		
		@inbounds for i = 2:nx+1, j = 2:ny+1, k = 2:nz+1
					divu = ( u[i,j,k]-u[i-1,j,k] )/dx[i,j]
					divv = ( v[i,j,k]-v[i,j-1,k] )/dy[i,j]
					divw = ( w[i,j,k]-w[i,j,k+1] )/dz[k]
					div  = dtin*B[i,j,k]*( divu+divv+divw )
					
					khe = 0.5*( ah[i,j,k]+ah[i+1,j,k] )
					dife = khe*( rho[i+1,j,k]-rho[i,j,k] )/dx[i,j]
					if ( dry[i+1,j,k] != 0 ) dife = 0.0 end
					khw = 0.5*( ah[i,j,k]+ah[i-1,j,k] )
					difw = khw*( rho[i,j,k]-rho[i-1,j,k] )/dx[i,j]
					if ( dry[i-1,j,k] != 0 ) difw = 0.0 end
					khn = 0.5*( ah[i,j,k]+ah[i,j+1,k] )
					difn = khn*( rho[i,j+1,k]-rho[i,j,k] )/dy[i,j]
					if ( dry[i,j+1,k] != 0 ) difn = 0.0 end
					khs = 0.5*( ah[i,j,k]+ah[i,j-1,k] )
					difs = khs*( rho[i,j,k]-rho[i,j-1,k] )/dy[i,j]
					if ( dry[i,j-1,k] != 0 ) difs = 0.0 end
					
					difhor = (dife-difw)/dx[i,j] + (difn-difs)/dy[i,j]
					
					aztop = 0.5*( kz[i,j,k]+kz[i,j,k-1] )
					dift = aztop*( rho[i,j,k-1]-rho[i,j,k] )/dz[k]
					if ( dry[i,j,k-1] != 0 ) dift = 0.0 end
					azbot = 0.5*( kz[i,j,k]+kz[i,j,k+1] )
					difb = azbot*( rho[i,j,k]-rho[i,j,k+1] )/dz[k]
					if ( dry[i,j,k+1] != 0 ) difb = 0.0 end
					difver = ( dift-difb) /dz[k]
					dif = dtin*( difhor+difver )
					rhon[i,j,k] = rho[i,j,k] + BN[i,j,k] + div + dif
		end
		
		@inbounds for k = 1:nz+2, j = 1:ny+2
				rhon[1,j,k] = rhon[2,j,k]
				rhon[nx+2,j,k] = rhon[nx+1,j,k]
		end
		
		@inbounds for k = 1:nz+2, i = 1:nx+2
				rhon[i,1,k] = rhon[i,2,k]
				rhon[i,ny+2,k] = rhon[i,ny+1,k]
		end
		
		@inbounds for j = 1:ny+2, i = 1:nx+2
				rhon[i,j,1] = rhon[i,j,2]
				rhon[i,j,nz+2] = rhon[i,j,nz+1]
		end
		
		@. rho = rhon
		
    end