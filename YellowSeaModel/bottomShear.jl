# bottom shear stress

	function BottomShear(nx,ny,ib,drag,u,v,taubu,taubv)
		
		umean = 0.0
		vmean = 0.0
		speed = 0.0
		
		@inbounds for i = 2:nx+1, j = 2:ny+1
				umean = u[i,j,ib[i,j]]
				vmean = 0.25*( v[i,j,ib[i,j]] + v[i+1,j,ib[i,j]] + v[i,j-1,ib[i,j]] + v[i+1,j-1,ib[i,j]] )
				speed = sqrt( umean*umean + vmean*vmean )
				taubu[i,j] = drag * umean * speed
				
				vmean = v[i,j,ib[i,j]]
				umean = 0.25*( u[i,j,ib[i,j]] + u[i-1,j,ib[i,j]] + u[i,j+1,ib[i,j]] + u[i-1,j+1,ib[i,j]] )
				speed = sqrt( umean*umean + vmean*vmean )
				taubv[i,j] = drag * vmean * speed
		end
		
    end