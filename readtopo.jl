# read topo.csv

	function ReadTopo(nx,ny,nz,depth,ele,ib,dx,dy,land,dry)

		fname = "./input/topo.csv"
		rdata = readcsv(fname)
		l=0
		@inbounds for j = 1:ny+2, i = 1:nx+2
				l = l+1
				dx[i,j] = rdata[l,3]
				dy[i,j] = rdata[l,4]
				ele[i,j] = rdata[l,5]
				depth[i,j] = max(-ele[i,j],0)
					if (dx[i,j] == 0.0)||(dy[i,j] == 0.0) 
 						land[i,j] = 1
 						dry[i,j,2:nz+2] = 1
 					end
 				na = trunc(Int64,depth[i,j])
 				for k = 2:nz+1
					nb = k-1
 					zsum = sum(dz[2:k])
 					if ( na <= zsum )
	 					@goto finish
 					end
	 			end
	 			#
	 			@label finish
	 			nb = min(nb,nz)
 				for k = nb+2:nz+2	# this point is very important nb+2, nz+2
					dry[i,j,k] = 1
				end
		end
			
		#fname = "./input/delxy.csv"
		#rdata = readcsv(fname)
		#l=0
		@inbounds for j = 1:ny+2, i = 1:nx+2
		#		l = l+1
		#		dx[i,j] = rdata[l,5]
		#		dy[i,j] = rdata[l,6]
				# reset the ib(bottom index)
				for k = 1:nz+1
					if (dry[i,j,k] == 0) && (dry[i,j,k+1] != 0) ib[i,j] = k end
				end
		end
				
	end