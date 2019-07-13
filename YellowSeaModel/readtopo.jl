# read topo.csv

	function ReadTopo(depth,ele,ib,dx,dy,land,dry)

		fname = "./input/topo.csv"
		rdata = readcsv(fname)
		l=0
		@inbounds for j = 1:ny+2, i = 1:nx+2
				l = l+1
				depth[i,j] = -rdata[l,6]
				ele[i,j] = rdata[l,6]
					if (depth[i,j] <= 0.0)
 						land[i,j] = 1
 						dry[i,j,2:nz+2] = 1
 					end
				zsum = 0.0
 				na = trunc(Int64,depth[i,j])
 				for k = 2:nz+1
 					nb = k-1
	 				zsum = zsum+dz[k]
	 				if ( land[i,j]==0 )
		 				dry[i,j,k] = 0.0
		 			end
	 				if ( na <= zsum )
		 				@goto finish
	 				end
	 			end
	 				
	 			@label finish
 				nb = min(nb,nz)
 				for k = nb+2:nz+2	# this point is very important nb+2, nz+2
					dry[i,j,k] = 1
				end
				
				# reset the ib(bottom index)
				for k = 1:nz+1
					if (dry[i,j,k] == 0) && (dry[i,j,k+1] != 0) ib[i,j] = k end
				end
		end
		
	end