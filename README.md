# ePnP_eigen


ePnP算法的eigen实现

不BB，三行代码即可调用：

1：生成epnp类 PnP；
2：输入设定参数：PnP.PnP.set_parameter(n, P3, P2, iM);
3：输出计算R和t：PnP.compute_pose(R_est, t_est)



还有：PnP.print_pose(R_est, t_est);打印出来
float err2 = PnP.compute_pose(R_est, t_est)；得到投影误差。

还有：经过检验速度再release版本为0.15ms。

还有：设定参数的储存顺序说明：
n is point number	
point_3d storage as XW1,YW1,ZW1,XW2,YW2,ZW2.....Size 3*n
point_2d storage as u1, v1, u2, v2.....Size 2*n
inter_matrix storage as fx,fy,u0,v0 Size 4*1

R_est:double类型[3][3]数组，表示旋转
t_est:double类型[3]向量，表示位移

期望后续


