# advection

	function Advection(nx,ny,nz,B,BN,CuPos,CuNeg,CvPos,CvNeg,CwPos,CwNeg)
		
		RxPos = zeros(Float64, nx+2,ny+2,nz+2)	# ( B(i,j,k)-B(i-1,j,k) )/delB
		RxNeg = zeros(Float64, nx+2,ny+2,nz+2)	# ( B(i,j,k)-B(i,j-1,k) )/delB
		RyPos = zeros(Float64, nx+2,ny+2,nz+2)	# ( B(i,j,k)-B(i,j,k+1) )/delB
		RyNeg = zeros(Float64, nx+2,ny+2,nz+2)	# ( B(i+2,j,k)-B(i+1,j,k) )/delB
		RzPos = zeros(Float64, nx+2,ny+2,nz+2)	# ( B(i,j+2,k)-B(i,j+1,k) )/delB
		RzNeg = zeros(Float64, nx+2,ny+2,nz+2)	# ( B(i,j,k-2)-B(i,j,k-1) )/delB
	
		delB = 0.0
		term1 = 0.0
		term2 = 0.0
		term3 = 0.0
		term4 = 0.0
		term5 = 0.0
		term6 = 0.0
		BwPos = 0.0
		BwNeg = 0.0
		BePos = 0.0
		BeNeg = 0.0
		BnPos = 0.0
		BnNeg = 0.0
		BsPos = 0.0
		BsNeg = 0.0
		BbPos = 0.0
		BbNeg = 0.0
		BtPos = 0.0
		BtNeg = 0.0
		R = 0.0
		
		@fastmath @inbounds for k = 2:nz+1, j = 2:ny+1, i = 2:nx+1
					delB = B[i+1,j,k] - B[i,j,k]
					if ( abs(delB)>0.0 ) RxPos[i,j,k] = ( B[i,j,k]-B[i-1,j,k] )/delB end
					delB = B[i,j+1,k] - B[i,j,k]
					if ( abs(delB)>0.0 ) RyPos[i,j,k] = ( B[i,j,k]-B[i,j-1,k] )/delB end
					delB = B[i,j,k-1] - B[i,j,k]
					if ( abs(delB)>0.0 ) RzPos[i,j,k] = ( B[i,j,k]-B[i,j,k+1] )/delB end
		end
		
		@fastmath @inbounds for k = 2:nz+1, j = 2:ny+1, i = 1:nx
					delB = B[i+1,j,k] - B[i,j,k]
					if ( abs(delB)>0.0 ) RxNeg[i,j,k] = ( B[i+2,j,k]-B[i+1,j,k] )/delB end
		end
		
		@fastmath @inbounds for k = 2:nz+1, j = 1:ny, i = 2:nx+1
					delB = B[i,j+1,k] - B[i,j,k]
					if ( abs(delB)>0.0 ) RyNeg[i,j,k] = ( B[i,j+2,k]-B[i,j+1,k] )/delB end
		end
		
		@fastmath @inbounds for k = 3:nz+2, j = 2:ny+1, i = 2:nx+1
					delB = B[i,j,k-1] - B[i,j,k]
					if ( abs(delB)>0.0 ) RzNeg[i,j,k] = ( B[i,j,k-2]-B[i,j,k-1] )/delB end
		end
		
		@fastmath @inbounds for k = 2:nz+1, j = 2:ny+1, i = 2:nx+1
				
					term1 = (1.0-CuPos[i-1,j,k]) * (B[i,j,k]-B[i-1,j,k])
					R = RxPos[i-1,j,k]
					BwPos = B[i-1,j,k] + 0.5*Psi(R) * term1
					
					term1 = (1.0+CuNeg[i-1,j,k]) * (B[i,j,k]-B[i-1,j,k])
					R = RxNeg[i-1,j,k]
					BwNeg = B[i,j,k] - 0.5*Psi(R) * term1
					
					term1 = (1.0-CuPos[i,j,k]) * (B[i+1,j,k]-B[i,j,k])
					R = RxPos[i,j,k]
					BePos = B[i,j,k] + 0.5*Psi(R) * term1
					
					term1 = (1.0+CuNeg[i,j,k]) * (B[i+1,j,k]-B[i,j,k])
					R = RxNeg[i,j,k]
					BeNeg = B[i+1,j,k] - 0.5*Psi(R) * term1
					
					term1 = (1.0-CvPos[i,j-1,k]) * (B[i,j,k]-B[i,j-1,k])
					R = RyPos[i,j-1,k]
					BsPos = B[i,j-1,k] + 0.5*Psi(R) * term1
					
					term1 = (1.0+CvNeg[i,j-1,k]) * (B[i,j,k]-B[i,j-1,k])
					R = RyNeg[i,j-1,k]
					BsNeg = B[i,j,k] - 0.5*Psi(R) * term1
					
					term1 = (1.0-CvPos[i,j,k]) * (B[i,j+1,k]-B[i,j,k])
					R = RyPos[i,j,k]
					BnPos = B[i,j,k] + 0.5*Psi(R) * term1
					
					term1 = (1.0+CvNeg[i,j,k]) * (B[i,j+1,k]-B[i,j,k])
					R = RyNeg[i,j,k]
					BnNeg = B[i,j+1,k] - 0.5*Psi(R) * term1
					
					term1 = (1.0-CwPos[i,j,k+1]) * (B[i,j,k]-B[i,j,k+1])
					R = RzPos[i,j,k+1]
					BbPos = B[i,j,k+1] + 0.5*Psi(R) * term1
					
					term1 = (1.0+CwNeg[i,j,k+1]) * (B[i,j,k]-B[i,j,k+1])
					R = RzNeg[i,j,k+1]
					BbNeg = B[i,j,k] - 0.5*Psi(R) * term1
					
					term1 = (1.0-CwPos[i,j,k]) * (B[i,j,k-1]-B[i,j,k])
					R = RzPos[i,j,k]
					BtPos = B[i,j,k] + 0.5*Psi(R) * term1
					
					term1 = (1.0+CwNeg[i,j,k]) * (B[i,j,k-1]-B[i,j,k])
					R = RzNeg[i,j,k]
					BtNeg = B[i,j,k-1] - 0.5*Psi(R) * term1
					
					term1 = (CuPos[i-1,j,k]*BwPos) + (CuNeg[i-1,j,k]*BwNeg)
					term2 = (CuPos[i,j,k]*BePos) + (CuNeg[i,j,k]*BeNeg)
					term3 = (CvPos[i,j-1,k]*BsPos) + (CvNeg[i,j-1,k]*BsNeg)
					term4 = (CvPos[i,j,k]*BnPos) + (CvNeg[i,j,k]*BnNeg)
					term5 = (CwPos[i,j,k+1]*BbPos) + (CwNeg[i,j,k+1]*BbNeg)
					term6 = (CwPos[i,j,k]*BtPos) + (CwNeg[i,j,k]*BtNeg)
					
					BN[i,j,k] = term1-term2 + term3-term4 + term5-term6
		end
		
    end
    
    
    function Psi(R)
    
    	R = max( min(2.0*R ,1.0), min(R, 2.0), 0.0)
    	
    end