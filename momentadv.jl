# non-linear term. momentum advection.

	function MomentumAdvection(nx,ny,nz,dx,dy,dz,dt,dtin,u,v,w,B,BN,advx,advy,advz)
		#@sync begin
			AdU(nx,ny,nz,dx,dy,dz,dt,dtin,u,v,w,B,BN,advx)
			AdV(nx,ny,nz,dx,dy,dz,dt,dtin,u,v,w,B,BN,advy)
			AdW(nx,ny,nz,dx,dy,dz,dt,dtin,u,v,w,B,BN,advz)
		#end
		
	end
	
	function AdU(nx,ny,nz,dx,dy,dz,dt,dtin,u,v,w,B,BN,advx)
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
		
		ex1 = 0.0
		ex2 = 0.0
		ex3 = 0.0
		ex4 = 0.0
		ex5 = 0.0
		ex6 = 0.0
		
		# u-momentum
		@inbounds for k = 1:nz+2, j = 1:ny+1, i = 1:nx+1
					ex1 = u[i,j,k]+u[i+1,j,k] + abs(u[i,j,k]) + abs(u[i+1,j,k])
					ex2 = u[i,j,k]+u[i+1,j,k] - abs(u[i,j,k]) - abs(u[i+1,j,k])
					ex3 = v[i,j,k]+v[i+1,j,k] + abs(v[i,j,k]) + abs(v[i+1,j,k])
					ex4 = v[i,j,k]+v[i+1,j,k] - abs(v[i,j,k]) - abs(v[i+1,j,k])
					ex5 = w[i,j,k]+w[i+1,j,k] + abs(w[i,j,k]) + abs(w[i+1,j,k])
					ex6 = w[i,j,k]+w[i+1,j,k] - abs(w[i,j,k]) - abs(w[i+1,j,k])
					CuPos[i,j,k] = 0.25*( ex1 )*dtin/dx[i,j]
					CuNeg[i,j,k] = 0.25*( ex2 )*dtin/dx[i,j]
					CvPos[i,j,k] = 0.25*( ex3 )*dtin/dy[i,j]
					CvNeg[i,j,k] = 0.25*( ex4 )*dtin/dy[i,j]
					CwPos[i,j,k] = 0.25*( ex5 )*dtin/dz[k]
					CwNeg[i,j,k] = 0.25*( ex6 )*dtin/dz[k]
		end
		
		@. B = u
		
		Advection(nx,ny,nz,B,BN,CuPos,CuNeg,CvPos,CvNeg,CwPos,CwNeg)
		
		@inbounds for k = 2:nz+1, j = 2:ny+1, i = 2:nx+1
					ex1 = u[i+1,j,k] - u[i-1,j,k]
					ex2 = v[i,j,k] + v[i+1,j,k] - v[i,j-1,k] - v[i+1,j-1,k]
					ex3 = w[i,j,k] + w[i+1,j,k] - w[i,j,k+1] - w[i+1,j,k+1]
					divu = 0.5*( ex1 )/dx[i,j]
					divv = 0.5*( ex2 )/dy[i,j]
					divw = 0.5*( ex3 )/dz[k]
					div = dtin*B[i,j,k]*( divu+divv+divw )
					advx[i,j,k] = BN[i,j,k] + div
		end
		
	end
		
		
	function AdV(nx,ny,nz,dx,dy,dz,dt,dtin,u,v,w,B,BN,advy)
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
		
		ex1 = 0.0
		ex2 = 0.0
		ex3 = 0.0
		ex4 = 0.0
		ex5 = 0.0
		ex6 = 0.0
		
		# v-momentum
		@inbounds for k = 1:nz+2, j = 1:ny+1, i = 1:nx+1
					ex1 = u[i,j,k]+u[i,j+1,k] + abs(u[i,j,k]) + abs(u[i,j+1,k])
					ex2 = u[i,j,k]+u[i,j+1,k] - abs(u[i,j,k]) - abs(u[i,j+1,k])
					ex3 = v[i,j,k]+v[i,j+1,k] + abs(v[i,j,k]) + abs(v[i,j+1,k])
					ex4 = v[i,j,k]+v[i,j+1,k] - abs(v[i,j,k]) - abs(v[i,j+1,k])
					ex5 = w[i,j,k]+w[i,j+1,k] + abs(w[i,j,k]) + abs(w[i,j+1,k])
					ex6 = w[i,j,k]+w[i,j+1,k] - abs(w[i,j,k]) - abs(w[i,j+1,k])
					CuPos[i,j,k] = 0.25*( ex1 )*dtin/dx[i,j]
					CuNeg[i,j,k] = 0.25*( ex2 )*dtin/dx[i,j]
					CvPos[i,j,k] = 0.25*( ex3 )*dtin/dy[i,j]
					CvNeg[i,j,k] = 0.25*( ex4 )*dtin/dy[i,j]
					CwPos[i,j,k] = 0.25*( ex5 )*dtin/dz[k]
					CwNeg[i,j,k] = 0.25*( ex6 )*dtin/dz[k]
		end
		
		@. B = v
		
		Advection(nx,ny,nz,B,BN,CuPos,CuNeg,CvPos,CvNeg,CwPos,CwNeg)
		
		@inbounds for k = 2:nz+1, j = 2:ny+1, i = 2:nx+1
					ex1 = u[i,j,k] + u[i,j+1,k] - u[i-1,j,k] - u[i-1,j+1,k]
					ex2 = v[i,j+1,k] - v[i,j-1,k]
					ex3 = w[i,j,k] + w[i,j+1,k] - w[i,j,k+1] - w[i,j+1,k+1]
					divu = 0.5*( ex1 )/dx[i,j]
					divv = 0.5*( ex2 )/dy[i,j]
					divw = 0.5*( ex3 )/dz[k]
					div = dtin*B[i,j,k]*( divu+divv+divw )
					advy[i,j,k] = BN[i,j,k] + div
		end
	end
	
	
	function AdW(nx,ny,nz,dx,dy,dz,dt,dtin,u,v,w,B,BN,advz)
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
		
		ex1 = 0.0
		ex2 = 0.0
		ex3 = 0.0
		ex4 = 0.0
		ex5 = 0.0
		ex6 = 0.0
		
		# w-momentum
		@inbounds for k = 1:nz+1, j = 1:ny+2, i = 1:nx+2
						ex1 = u[i,j,k]+u[i,j,k-1] + abs(u[i,j,k]) + abs(u[i,j,k-1])
						ex2 = u[i,j,k]+u[i,j,k-1] - abs(u[i,j,k]) - abs(u[i,j,k-1])
						ex3 = v[i,j,k]+v[i,j,k-1] + abs(v[i,j,k]) + abs(v[i,j,k-1])
						ex4 = v[i,j,k]+v[i,j,k-1] - abs(v[i,j,k]) - abs(v[i,j,k-1])
						ex5 = w[i,j,k]+w[i,j,k-1] + abs(w[i,j,k]) + abs(w[i,j,k-1])
						ex6 = w[i,j,k]+w[i,j,k-1] - abs(w[i,j,k]) - abs(w[i,j,k-1])
						CuPos[i,j,k] = 0.25*( ex1 )*dtin/dx[i,j]
						CuNeg[i,j,k] = 0.25*( ex2 )*dtin/dx[i,j]
						CvPos[i,j,k] = 0.25*( ex3 )*dtin/dy[i,j]
						CvNeg[i,j,k] = 0.25*( ex4 )*dtin/dy[i,j]
						CwPos[i,j,k] = 0.25*( ex5 )*dtin/dz[k]
						CwNeg[i,j,k] = 0.25*( ex6 )*dtin/dz[k]
		end
		
		@. B = w
		
		Advection(nx,ny,nz,B,BN,CuPos,CuNeg,CvPos,CvNeg,CwPos,CwNeg)
		
		@inbounds for k = 2:nz+1, j = 2:ny+1, i = 2:nx+1
					ex1 = u[i,j,k] + u[i,j,k-1] - u[i-1,j,k] - u[i-1,j,k-1]
					ex2 = v[i,j,k] + v[i,j,k-1] - v[i,j-1,k] - v[i,j-1,k-1]
					ex3 = w[i,j,k-1] - w[i,j,k+1]
					divu = 0.5*( ex1 )/dx[i,j]
					divv = 0.5*( ex2 )/dy[i,j]
					divw = 0.5*( ex3 )/dz[k]
					div = dtin*B[i,j,k]*( divu+divv+divw )
					advz[i,j,k] = BN[i,j,k] + div
		end
		
    end