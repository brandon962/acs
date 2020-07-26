github : https://github.com/brandon962/acs

1. 如何執行

		a. 
		執行 start.bat 
		將預設參數並畫出最短路線圖及收斂曲線圖
		預設參數設定為 :
			runs       : 30
			iterations : 1000
			nodes      : 51
			ants       : 51
			pheromone  : 1
			alpha      : 1
			beta       : 2.8
			rho        : 0.9
			file input : eil51.txt

		b.
		[main.exe] [runs] [iterations] [nodes] [ants] [pheromone] [alpha] [beta] [rho] [file_input]

2. 執行結果
	
		執行30 run 
		每run平均最短路徑 : 432.807
		曾發現最短路徑 : 428.872
		總花費時間 : 34 秒

![image](https://github.com/brandon962/acs/blob/master/path.png)

3. 收斂圖

![image](https://github.com/brandon962/acs/blob/master/convergence.png)


