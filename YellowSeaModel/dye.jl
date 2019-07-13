# density advection

	function Dye(nx,ny,nz,dry,dx,dy,dz,dt,dtin,dye,dyen,B,BN,u,v,w,ah,kz,ib,CuPos,CuNeg,CvPos,CvNeg,CwPos,CwNeg)
		
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
					CuPos[i,j,k] = 0.5*( u[i,j,k] + abs(u[i,j,k]) ) * dtin[1]/dx[i,j]
					CuNeg[i,j,k] = 0.5*( u[i,j,k] - abs(u[i,j,k]) ) * dtin[1]/dx[i,j]
					CvPos[i,j,k] = 0.5*( v[i,j,k] + abs(v[i,j,k]) ) * dtin[1]/dy[i,j]
					CvNeg[i,j,k] = 0.5*( v[i,j,k] - abs(v[i,j,k]) ) * dtin[1]/dy[i,j]
					CwPos[i,j,k] = 0.5*( w[i,j,k] + abs(w[i,j,k]) ) * dtin[1]/dz[k]
					CwNeg[i,j,k] = 0.5*( w[i,j,k] - abs(w[i,j,k]) ) * dtin[1]/dz[k]
		end
		
		@. B = dye
		
		Advection(nx,ny,nz,B,BN,CuPos,CuNeg,CvPos,CvNeg,CwPos,CwNeg)
		
		@inbounds for i = 2:nx+1, j = 2:ny+1, k = 2:nz+1
				if (dry[i,j,k] == 0 )
					divu = ( u[i,j,k]-u[i-1,j,k] )/dx[i,j]
					divv = ( v[i,j,k]-v[i,j-1,k] )/dy[i,j]
					divw = ( w[i,j,k]-w[i,j,k+1] )/dz[k]
					div  = dtin[1]*B[i,j,k]*( divu+divv+divw )
					
					khe = 0.5*( ah[i,j,k]+ah[i+1,j,k] )
					dife = khe*( dye[i+1,j,k]-dye[i,j,k] )/dx[i,j]
					if ( dry[i+1,j,k] != 0 ) dife = 0.0 end
					khw = 0.5*( ah[i,j,k]+ah[i-1,j,k] )
					difw = khw*( dye[i,j,k]-dye[i-1,j,k] )/dx[i,j]
					if ( dry[i-1,j,k] != 0 ) difw = 0.0 end
					khn = 0.5*( ah[i,j,k]+ah[i,j+1,k] )
					difn = khn*( dye[i,j+1,k]-dye[i,j,k] )/dy[i,j]
					if ( dry[i,j+1,k] != 0 ) difn = 0.0 end
					khs = 0.5*( ah[i,j,k]+ah[i,j-1,k] )
					difs = khs*( dye[i,j,k]-dye[i,j-1,k] )/dy[i,j]
					if ( dry[i,j-1,k] != 0 ) difs = 0.0 end
					
					difhor = (dife-difw)/dx[i,j] + (difn-difs)/dy[i,j]
					
					aztop = 0.5*( kz[i,j,k]+kz[i,j,k-1] )
					dift = aztop*( dye[i,j,k-1]-dye[i,j,k] )/dz[k]
					if ( dry[i,j,k-1] != 0 ) dift = 0.0 end
					azbot = 0.5*( kz[i,j,k]+kz[i,j,k+1] )
					difb = azbot*( dye[i,j,k]-dye[i,j,k+1] )/dz[k]
					if ( dry[i,j,k+1] != 0 ) difb = 0.0 end
					difver = ( dift-difb) /dz[k]
					dif = dtin[1]*( difhor+difver )
					dyen[i,j,k] = dye[i,j,k] + BN[i,j,k] + div + dif
				end
		end
		
		@inbounds for k = 1:nz+2, j = 1:ny+2
				dyen[1,j,k] = dyen[2,j,k]
				dyen[nx+2,j,k] = dyen[nx+1,j,k]
		end
		
		@inbounds for k = 1:nz+2, i = 1:nx+2
				dyen[i,1,k] = dyen[i,2,k]
				dyen[i,ny+2,k] = dyen[i,ny+1,k]
		end
		
		@inbounds for j = 1:ny+2, i = 1:nx+2
				dyen[i,j,1] = dyen[i,j,2]
				dyen[i,j,nz+2] = dyen[i,j,nz+1]
		end
		
		@. dye = dyen
		
    end