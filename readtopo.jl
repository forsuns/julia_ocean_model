# read topo.csv

	function ReadTopo(ele,ib)

		fname = "./input/topo.csv"
		rdata = readcsv(fname)
		k=0
		@inbounds for j = 1:ny+2, i = 1:nx+2
				k = k+1
				ele[i,j] = rdata[k,3]
				ib[i,j] = 1
		end
		
	end