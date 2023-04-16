# two-stage-robust-MG

复现中国电机工程学报《微电网两阶段鲁棒优化经济调度方法》

根据文中的强对偶理论编程求解时出现了一些问题，因此重新推导了模型的KKT条件进行求解

语言：Python 3.10.1 + Gurobi 10.0.1

程序说明：twostageMG.py为非紧凑形式的约束，KKTmatrix.py将非紧凑形式的约束转化为紧凑形式，MGCCGKKT为采用KKT方法的CCG两阶段鲁棒求解程序，运行MGCCGKKT.py即可

![image](https://user-images.githubusercontent.com/51228607/232202948-6b38c3f2-0d30-403a-bfc7-0d0c1014b106.png)

可控分布式电源出力

![image](https://user-images.githubusercontent.com/51228607/232203153-e5c4c9cf-462e-41c5-a436-24bc2c8ca893.png)

光伏出力

![image](https://user-images.githubusercontent.com/51228607/232203168-4ec9a185-1041-4f85-b580-3993cd200758.png)

负荷
