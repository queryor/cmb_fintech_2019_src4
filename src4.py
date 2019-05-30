# 任务1
import numpy as np
x0 = np.array([22.5369548317,114.0225501220])
x1 = np.array([39.9064024197,116.3624819407])
x2 = np.array([36.0541745829,103.8366184893])
x3 = np.array([31.2427704453,121.5058705816])
x4 = np.array([45.7680708990, 126.6203889331])

from geographiclib.geodesic import Geodesic
a = Geodesic(6378137,1/298.257222101)
i1 = a.WGS84.Inverse(x1[0],x1[1],x0[0],x0[1])['s12']
i2 = a.WGS84.Inverse(x2[0],x2[1],x0[0],x0[1])['s12']
i3 = a.WGS84.Inverse(x3[0],x3[1],x0[0],x0[1])['s12']
i4 = a.WGS84.Inverse(x4[0],x4[1],x0[0],x0[1])['s12']
i1,i2,i3,i4 = round(i1,5),round(i2,5),round(i3,5),round(i4,5)
print(i1,i2,i3,i4)
# 任务一计算结果
data1 = [i1,i2,i3,i4]

# 任务2 Dijkstra算法
# 任务二、任务三的数据读取
import geopandas
import matplotlib.pyplot as plt
import numpy as np
%matplotlib inline

n_path = 'data/FT_Camp_4/sandiego/nodes/nodes.shp' 
node_shp_df = geopandas.GeoDataFrame.from_file(n_path) #读取nodes.shp数据
e_path = 'data/FT_Camp_4/sandiego/edges/edges.shp' 
edge_shp_df = geopandas.GeoDataFrame.from_file(e_path) #读取edges.shp数据

begin=[]
begin.append([32.7496294, -117.1092834])
begin.append([32.740292, -117.103363])
begin.append([32.795853, -117.174637])
begin.append([32.740292, -117.103363])
begin.append([32.733253, -117.054108])
begin.append([32.760967, -117.02327])
begin.append([32.628593, -117.016388])
begin.append([32.7507351, -117.1092739])
begin.append([32.766434, -117.081886])
begin.append([32.681469, -116.992996])
end = []
end.append([32.755474, -117.096519])
end.append([32.789948, -117.240753])
end.append([32.623028, -117.070381])
end.append([32.73925, -117.133629])
end.append([32.733768, -116.965523])
end.append([32.731846, -117.032837])
end.append([32.6183834, -117.1340124])
end.append([32.612331, -117.06707])
end.append([32.740261, -117.096588])
end.append([32.790749, -116.9315478])



node_num = len(node_shp_df.osmid)
print("图中共有{}个点\n".format(node_num))
edge_num = len(edge_shp_df.osmid)
print("图中共有{}条边\n".format(edge_num))
Node_crd = node_shp_df["geometry"]

# 节点字典
node_dict = {}
Node = node_shp_df.osmid
for i,nid in enumerate(node_shp_df.osmid):
    node_dict[nid]=i
# 权重矩阵
## 边没有连接weight 为 inf
max_float=float('inf')
W = {}
print("W 初始化成功")

## 自己到自己的weight 为 0
# for i in range(node_num):
#     W[i][i]=0


## 记录每条边的索引 可以根据两个端点找到相应的边 如果是-1说明没有边
E = {}
print("E 初始化成功")
#记录每个点的邻近点
Node_near = [[] for i in range(node_num)]

from_list = [node_dict[i] for i in edge_shp_df["from"]]
to_list = [node_dict[i] for i in edge_shp_df["to"]]
oneway_list = edge_shp_df["oneway"]
length_list = [float(i) for i in edge_shp_df["length"]]
speed_list = [float(i) for i in edge_shp_df["base_speed"]]
weight_list = np.array(length_list)/np.array(speed_list)
osmid_list = edge_shp_df["osmid"]

## 根据每条边初始化相应的weight
for i in range(edge_num):
    W[str(from_list[i])+'_'+str(to_list[i])]=weight_list[i]
    E[str(from_list[i])+'_'+str(to_list[i])]=i
    Node_near[from_list[i]].append(to_list[i])
    if oneway_list[i] == "False":
        W[str(to_list[i])+'_'+str(from_list[i])]=weight_list[i]
        E[str(to_list[i])+'_'+str(from_list[i])]=i
        Node_near[to_list[i]].append(from_list[i])

max_float=float('inf')
## dijkstra算法实现
def Dijkstra(points,graph,start,end):
   # map = [[ max_float for i in range(points)] for j in range(points)]
    pre = [0]*(points) #记录前驱
    vis = [0]*(points) #记录节点遍历状态
    dis = [max_float for i in range(points)] #保存最短距离
    road = [0]*(points) #保存最短路径
    roads = []
    map = graph
 
    for i in range(points):#初始化起点到其他点的距离
        if i == start :
            dis[i] = 0
        else:
            key = str(start)+'_'+str(i)
            if key in map:
                dis[i]=map[key]
            else:
                dis[i] = max_float
        if dis[i] != max_float:
            pre[i] = start
        else:
            pre[i]=-1
    pre[start] = -1
    vis[start] = 1
    for i in range(points):#每循环一次确定一条最短路
        min = max_float
        for j in range(points):#寻找当前最短路
            if vis[j] == 0 and dis[j] < min :
                t = j
                min = dis[j]
        vis[t] = 1 #找到最短的一条路径 ,标记
        if t == end:
            break
        for j in Node_near[t]:
            key = str(t)+'_'+str(j)
            if key in map:
                if vis[j] == 0 and dis[j] > dis[t]+ map[key]:
                    dis[j] = dis[t] + map[key]
                    pre[j] = t
        
    p = end
    len = 0
    while p >= 0 and len < points:
        road[len] = p
        p = pre[p]
        len += 1
    mark = 0
    len -= 1
    while len >= 0:
        roads.append(road[len])
        len -= 1
    return dis[end],roads

## 起始id 和 目的id
## 得到起始id 和目的id
begin_id =[]
for crd in begin:
    for i,des in enumerate(Node_crd):
        if des.x == crd[1] and des.y ==crd[0]:
            begin_id.append(i)
print(begin_id)
end_id = []
for crd in end:
    for i,des in enumerate(Node_crd):
        if des.x == crd[1] and des.y ==crd[0]:
            end_id.append(i)
print(end_id)


## 计算最短路
dist =[]
roads = []
for i in range(len(begin_id)):
    dis,road = Dijkstra(node_num,W,int(begin_id[i]),int(end_id[i]))
    print(dis,road)
    dist.append(dis)
    roads.append(road)

# 任务二计算结果
data2 = []
for road in roads:
    t = []
    for i in road:
        t.append([Node_crd[i].y,Node_crd[i].x])
    data2.append(t)
#print(data2)

#任务三
# 1	(32.766434, -117.081886)
# 2	(32.7590688, -117.1618495)
# 3	(32.733253, -117.054108)
# 4	(32.623028, -117.070381)
CenterPoint=[[32.766434, -117.081886],[32.7590688, -117.1618495],[32.733253, -117.054108],[32.623028, -117.070381]]
## 得到center id
center_id =[]
for crd in CenterPoint:
    for i,des in enumerate(Node_crd):
        if des.x == crd[1] and des.y ==crd[0]:
            center_id.append(i)
print(center_id)

## 更新weight
W = {}
E = {}
for i in range(edge_num):
    W[str(from_list[i])+'_'+str(to_list[i])]=length_list[i]
    E[str(from_list[i])+'_'+str(to_list[i])]=i
    if oneway_list[i] == "False":
        W[str(to_list[i])+'_'+str(from_list[i])]=length_list[i]
        E[str(to_list[i])+'_'+str(from_list[i])]=i

## 给定图计算服务区域
import numpy as np
import matplotlib.pyplot as plt
import math


class point:
    def __init__(self,x,y):
        self.x = x
        self.y = y
def dis(a:point,b:point):
    return math.sqrt((a.x - b.x)**2 + (a.y - b.y)**2)

def get_node(Tree,V,id,pre_p,ans=[]):
    kids = Tree[id].copy()
    if len(kids)==0:
        ans.append(id)
        return
    pro = []
    for k in kids:
        dy = V[k].y-V[id].y
        dx = V[k].x-V[id].x
        p = math.atan2(dy,dx)
        p = (p-pre_p)%(2*math.pi)
        pro.append(p)
    pro = np.array(pro)
    kids = np.array(kids)
    index = pro.argsort()
    kids = kids[index]
    for k in kids:
        p = math.atan2(V[id].y-V[k].y,V[id].x-V[k].x)
        get_node(Tree,V,k,p,ans)
    

def isRayIntersectsSegment(poi,s_poi,e_poi): #[x,y] [lng,lat]
    #输入：判断点，边起点，边终点，都是[lng,lat]格式数组
    if s_poi[1]==e_poi[1]: #排除与射线平行、重合，线段首尾端点重合的情况
        return False
    if s_poi[1]>poi[1] and e_poi[1]>poi[1]: #线段在射线上边
        return False
    if s_poi[1]<poi[1] and e_poi[1]<poi[1]: #线段在射线下边
        return False
    if s_poi[1]==poi[1] and e_poi[1]>poi[1]: #交点为下端点，对应spoint
        return False
    if e_poi[1]==poi[1] and s_poi[1]>poi[1]: #交点为下端点，对应epoint
        return False
    if s_poi[0]<poi[0] and e_poi[0]<poi[0]: #线段在射线左边
        return False

    xseg=e_poi[0]-(e_poi[0]-s_poi[0])*(e_poi[1]-poi[1])/(e_poi[1]-s_poi[1]) #求交
    if xseg<poi[0]: #交点在射线起点的左侧
        return False
    return True  #排除上述情况之后

def isPoiWithinPoly(poi,poly):
    #输入：点，多边形二维数组
    #poly=[[x1,y1],[x2,y2],……,[xn,yn],[x1,y1]] 

    sinsc=0 #交点个数
    #循环每条边的曲线->each polygon 是二维数组[[x1,y1],…[xn,yn]]
    for i in range(len(poly)-1): #[0,len-1]
        s_poi=poly[i]
        e_poi=poly[i+1]
        if isRayIntersectsSegment(poi,s_poi,e_poi):
            sinsc+=1 #有交点就加1

    return True if sinsc%2==1 else  False

def dis_line(p1,pa,pb):
    x,y = p1.x,p1.y
    x1,y1,x2,y2 = pa.x,pa.y,pb.x,pb.y
    cross = (x2-x1)*(x-x1)+(y2-y1)*(y-y1)
    if cross<=0:
        return math.sqrt((x-x1)**2+(y-y1)**2)
    d2 = (x2-x1)**2 +(y2-y1)**2
    if cross>=d2:
        return math.sqrt((x-x2)**2+(y-y2)**2)
    r = cross/d2
    px = x1+(x2-x1)*r
    py = y1+(y2-y1)*r
    return math.sqrt((x-px)**2+(py-y)**2)

#凸包问题
import sys
import math
import time
import random
import datetime
#获取基准点的下标,基准点是p[k]
def get_leftbottompoint(p):
    k = 0
    for i in range(1, len(p)):
        if p[i][1] < p[k][1] or (p[i][1] == p[k][1] and p[i][0] < p[k][0]):
            k = i
    return k

#叉乘计算方法
def multiply(p1, p2, p0):
    return (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1])

#获取极角，通过求反正切得出，考虑pi/2的情况
def get_arc(p1, p0):
    # 兼容sort_points_tan的考虑
    if (p1[0] - p0[0]) == 0:
        if ((p1[1] - p0[1])) == 0:
            return -1
        else:
            return math.pi / 2
    tan = float((p1[1] - p0[1])) / float((p1[0] - p0[0]))
    arc = math.atan(tan)
    if (p1[1] - p0[1]) == 0 and float((p1[0] - p0[0]))<0:
        return math.pi
    if arc >= 0:
        return arc
    else:
        return math.pi + arc

#对极角进行排序,排序结果list不包含基准点
def sort_points_tan(p, pk):
    p2 = []
    for i in range(0, len(p)):
        p2.append({"index": i, "arc": get_arc(p[i], pk)})
    #print('排序前:',p2)
    p2.sort(key=lambda k: (k.get('arc')))
    #print('排序后:',p2)
    p_out = []
    for i in range(0, len(p2)):
        p_out.append(p[p2[i]["index"]])
    return p_out

def convex_hull(p):
    #p=list(set(p))
    #print('全部点:',p)
    k = get_leftbottompoint(p)
    pk = p[k]
    p.remove(p[k])
    #print('排序前去除基准点的所有点:',p,'基准点:',pk)
    p_sort = sort_points_tan(p, pk)   #按与基准点连线和x轴正向的夹角排序后的点坐标
    #print('其余点与基准点夹角排序:',p_sort)
    p_result = [pk,p_sort[0]]
    top = 2
    for i in range(1, len(p_sort)):
        #####################################
        #叉乘为正,向前递归删点;叉乘为负,序列追加新点
        while(multiply(p_result[-2], p_sort[i],p_result[-1]) > 0):
            p_result.pop()
        p_result.append(p_sort[i])    
    return p_result#测试


## 根据网络结果得到所有的路
def DFS_get_road(Tree,point,roads,road=[]):
    if len(Tree[point])==0:
        #print(road)
        roads.append(road.copy())
        return
    else:
        pro = []
        for j in Tree[point]:
            pro.append(get_arc([Node_crd[j].y,Node_crd[j].x],[Node_crd[point].y,Node_crd[point].x]))
        index = np.array(pro).argsort()
        #print(pro)
        nextpoints = np.array(Tree[point])[index].tolist().copy()
        for nextPoint in nextpoints:
            road.append(nextPoint)
            DFS_get_road(Tree,nextPoint,roads,road)
            road.pop()
            
# 测试
Tree = [[1,2,3],[4],[5],[6],[],[],[]]
all_roads = []
DFS_get_road(Tree,0,all_roads,[0])
print(all_roads)

def Dijkstra_dis(points,graph,distance,start,ans):
    ## points：点数
    ## garph: 权重矩阵
    ## start:起始点
    ## ans: 结果
    ## distance: 目标距离
   # map = [[ max_float for i in range(points)] for j in range(points)]
    pre = [0]*(points) #记录前驱
    vis = [0]*(points) #记录节点遍历状态
    dis = [max_float for i in range(points)] #保存最短距离
    map = graph
 
    for i in range(points):#初始化起点到其他点的距离
        if i == start :
            dis[i] = 0
        else :
            key = str(start)+'_'+str(i)
            if key in map:
                dis[i]=map[key]
            else:
                dis[i] = max_float
        if dis[i] != max_float:
            pre[i] = start
        else :
            pre[i] = -1
    vis[start] = 1
    for i in range(points):#每循环一次确定一条最短路
        min = max_float
        for j in range(points):#寻找当前最短路
            if vis[j] == 0 and dis[j] < min :
                t = j
                min = dis[j]
        vis[t] = 1 #找到最短的一条路径 ,标记
        if min>distance:
            break
        for j in Node_near[t]:
            key = str(t)+'_'+str(j)
            if key in map:
                if vis[j] == 0 and dis[j] > dis[t]+ map[key]:
                    dis[j] = dis[t] + map[key]
                    pre[j] = t
    all_point = []
    for i,l in enumerate(dis):
        if l>=distance and l<max_float:
            diffL = distance - dis[pre[i]]
            key = str(pre[i])+'_'+str(i)
            L = length_list[E[key]]
            pro = diffL/L
            crdX = (1-pro)*Node_crd[pre[i]].x + pro*Node_crd[i].x
            crdY = (1-pro)*Node_crd[pre[i]].y + pro*Node_crd[i].y
            ans.append([crdX,crdY])
            all_point.append([crdX,crdY])
               
    for i in range(node_num):
        if vis[i]==1:
            if dis[i]<distance:
                ans.append([Node_crd[i].x,Node_crd[i].y])
                all_point.append([Node_crd[i].x,Node_crd[i].y])
    ## 得到所有的线
    all_lines = []
    for i in range(node_num):
        if vis[i]==1:
            if dis[i]<distance and dis[i]>0:
                all_lines.append([(Node_crd[pre[i]].x,Node_crd[pre[i]].y),(Node_crd[i].x,Node_crd[i].y)])
                for j in Node_near[i]:
                    key = str(i)+'_'+str(j)
                    if length_list[E[key]]+dis[i]<=distance:
                        all_lines.append([(Node_crd[i].x,Node_crd[i].y),(Node_crd[j].x,Node_crd[j].y)])
                    else:
                        diffL = distance - dis[i]
                        key = str(i)+'_'+str(j)
                        L = length_list[E[key]]
                        pro = diffL/L
                        crdX = (1-pro)*Node_crd[i].x + pro*Node_crd[j].x
                        crdY = (1-pro)*Node_crd[i].y + pro*Node_crd[j].y
                        all_lines.append([(Node_crd[i].x,Node_crd[i].y),(crdX,crdY)])
    return all_point,all_lines
#         elif l<max_float and l>0:
#             crdX = Node_crd[i].x
#             crdY = Node_crd[i].y
#             ans.append([crdX,crdY])

# dj
data3 = []
all_point = []
all_roads = []
all_lines = []
for distance in [1000.0,3000.0]:
    for Point_id in center_id:
        print(Point_id)
        vis = [0 for i in range(node_num)]
        vis[Point_id]=1
        #DFS_dis(Point_id,vis,distance,ans)
        ans = []
        all_point_t,all_lines_t = Dijkstra_dis(node_num,W,distance,Point_id,ans)
        all_point.append(all_point_t)
        all_lines.append(all_lines_t)
        ans1 = convex_hull(ans)
        data3_t =[]
        for i in ans1:
            data3_t.append([i[1],i[0]])
        ## 加上起始点构成闭多边形
        data3_t.append([ans1[0][1],ans1[0][0]])
        data3.append(data3_t)
print("得到第三题答案")
# print(data3)

from shapely.geometry import LineString
from shapely.geometry import MultiLineString
lines = MultiLineString([LineString(i) for i in all_lines[1]])
ans = []
X = []
Y = []
data3 = []
for num in range(len(all_lines)):
    lines = MultiLineString([LineString(i) for i in all_lines[num]])
    b = lines.buffer(0.00093)
    try:
        bb = b.boundary[0]
    except TypeError:
        bb = b.boundary
    else:
        pass
    X_t = []
    Y_t = []
    ans_t = []
    for c in  bb.coords[:]:
        X_t.append(c[0])
        Y_t.append(c[1])
        ans_t.append([c[1],c[0]])
    X.append(X_t)
    Y.append(Y_t)
    data3.append(ans_t)

## 任务4
from sklearn.cluster import KMeans
from sklearn.externals import joblib
import numpy as np
import time
import matplotlib.pyplot as plt
import random
import math
# 任务四的数据读取
import pandas as pd
from tqdm import tqdm_notebook as tqdm
from concurrent.futures import ThreadPoolExecutor
df = pd.read_csv('data/FT_Camp_4/resident.csv')
df.head()
## step 1: 加载数据
print( "step 1: load data...")
dataSet_home = []
dataSet_company = []
for i in df['home']:
    x,y = i[1:-1].split(',')
    dataSet_home.append([float(y),float(x)])
for i in df['company']:
    x,y = i[1:-1].split(',')
    dataSet_company.append([float(y),float(x)])

dataSet_home = np.array(dataSet_home)
dataSet_company = np.array(dataSet_company)

# 固定半径的圆能覆盖的最多点数

class interval:
    def __init__(self,arg=0,flag=False):
        self.arg = arg
        self.flag = flag
class point:
    def __init__(self,x,y):
        self.x = x
        self.y = y
def dis(a:point,b:point):
    return math.sqrt((a.x - b.x)**2 + (a.y - b.y)**2)

# 画圆
def plot_circle(crd:point,r,axes,color='r'):
    #crd 圆心坐标
    # r  半径
    a,b = crd.x,crd.y
    theta = np.arange(0, 2*np.pi, 0.01)
    x = a + r * np.cos(theta)
    y = b + r * np.sin(theta)    
    axes.plot(x, y,color=color)
    axes.axis('equal')
def f1(p1,p2,event,r):
    dist = dis(p1,p2)
    if dist<=2.0*r:
        cta = math.atan2(p2.y-p1.y,p2.x-p1.x)
        if cta<0:
            cta+=2*math.pi
        delta = math.acos(dist/2/r)
        a1 = cta - delta
        a2 = cta + delta
        if a1<0:
            event.append(interval(arg=a1+2*math.pi,flag=True))
            event.append(interval(arg=a2+2*math.pi,flag=False))
        else:
            event.append(interval(arg=a1,flag=True))
            event.append(interval(arg=a2,flag=False))
            event.append(interval(arg=a1+2*math.pi,flag=True))
            event.append(interval(arg=a2+2*math.pi,flag=False))
def counter_point(points,center,r):
    counter = 0
    for p in points:
        if dis(center,p)<=r:
            counter+=1
    return counter
def  get_max_circle(points,r,name=""):
    d = datetime.datetime.now().strftime('%Y-%m-%d-%H-%M')
    f_name = "log/log_{}_{}.txt".format(name,d)
    
    n = len(points)
    # float的精度考虑
    ans_points = []
    distance={}
    #event = [interval() for i in range(4*n)]
    points.sort(key=lambda a:a.y)
    pbar = tqdm(total=n)
    start = time.time()
    for i in range(n):
        if i%10000==0 and i!=0:
            file = open(f_name,'a')
            ans_points.sort(key=lambda x:x[1],reverse=True)
            ans = []
            for k in range(len(ans_points)):
                p = ans_points[k][0]
                flag = 0
                for j in ans:
                    if dis(point(j[1],j[0]),p)<=0.04:
                        flag = 1
                        break
                if flag ==1:
                    continue
                else:
                    count = 0
                    for p1 in points:
                        if dis(p1,p)<=0.02:
                            count+=1
                    ans.append([p.y,p.x,count])
                    if len(ans)==5:
                        break
            ans.sort(key=lambda x:x[-1],reverse=True)
            file.write("{}/{}:{}    {}\n".format(i,n,time.time()-start,datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')))
            for s in ans:
                file.write(str(s)+'\n')
            file.close()
        pbar.update(1)
        event = []
        ans = -1
        #start = time.time()
        for j in range(n):
            if i == j:
                continue
            if j>i and dis(point(0,points[i].y),point(0,points[j].y))>2.0*r:
                break
            dist = dis(points[i],points[j])                
            if dist<=2.0*r:
                cta = math.atan2(points[j].y-points[i].y,points[j].x-points[i].x)
                if cta<0:
                    cta+=2*math.pi
                delta = math.acos(dist/2/r)
                a1 = cta - delta
                a2 = cta + delta
                if a1<0:
                    event.append(interval(arg=a1+2*math.pi,flag=True))
                    event.append(interval(arg=a2+2*math.pi,flag=False))
                else:
                    event.append(interval(arg=a1,flag=True))
                    event.append(interval(arg=a2,flag=False))
                    event.append(interval(arg=a1+2*math.pi,flag=True))
                    event.append(interval(arg=a2+2*math.pi,flag=False))
        if len(event)<ans:
            continue
        
        event.sort(key=lambda x:x.arg)
        #print("stage1:{}".format(time.time()-start))
        #print([x.flag for x in event_t[:num]]) 
        res = 0
        #print("stage2:{}".format(time.time()-start))
        #start = time.time()
        for j in range(len(event)):
            if event[j].flag:
                res+=1
            else: 
                res-=1
            if res>ans:
                arg = (event[j].arg+event[j+1 if j+1<len(event) else j].arg)/2
                p = points[i] 
                x = p.x + r*math.cos(arg)
                y = p.y + r*math.sin(arg)
                ans = res
                center_id = i
                center_arg = arg
        if ans!=-1:
            ans_points.append([point(x,y),ans])
        #print("stage3:{}".format(time.time()-start))
    p = points[center_id] 
    x = p.x + r*math.cos(center_arg)
    y = p.y + r*math.sin(center_arg)
    #print(len(ans_points))
    #print(ans_points)
    ans_points.sort(key=lambda x:x[1],reverse=True)
    #print(ans_points)
    return point(x,y),ans_points
    
## 测试
#画图初始化 
fig, axes = plt.subplots(figsize=(5,5))
# axes = fig.add_subplot(1,1,1)
##固定半径
r = 2
#初始化点坐标
x = [1,2,3,4,5,6,7]
y = [2,4,7,4,5,8,1]
axes.scatter(x,y)
points = [point(x[i],y[i]) for i in range(len(x))]
#画出原始圆
for p in points:
    plot_circle(p,r,axes)

## 得到覆盖点最多的圆的圆心
centerPoint,ans_points = get_max_circle(points,r)
plot_circle(ans_points[0][0],r,axes,color='b')
#显示图
#plt.show()      

# 计算公司
r = 0.02
# for p in AllPoints[:2000]:
#     plot_circle(p,r,axes)
AllPoints = [point(i[0],i[1]) for i in dataSet_company]
Tpoints = AllPoints
ans_company = []
i=0
random.shuffle(Tpoints)
print(len(Tpoints))
start = time.time()
centerPoint,ans_points = get_max_circle(Tpoints[:300000],r,name="comany")
print("time:{}".format(time.time()-start))
print(len(ans_points))
## 将结果 解析出来
ans_company = []
## 将结果 解析出来
pbar = tqdm(total=5)
for i in range(len(ans_points)):
    p = ans_points[i][0]
    flag = 0
    for j in ans_company:
        if dis(point(j[1],j[0]),p)<0.04:
            flag = 1
            break
    if flag ==1:
        continue
    else:
        count = 0
        for p1 in Tpoints:
            if dis(p1,p)<=0.02:
                count+=1
        ans_company.append([p.y,p.x,count])
        pbar.update(1)
        if len(ans_company)==5:
            break
ans_company.sort(key=lambda x:x[-1],reverse=True)
print(ans_company)


## 计算家庭
r = 0.02
# for p in AllPoints[:2000]:
#     plot_circle(p,r,axes)
AllPoints = [point(i[0],i[1]) for i in dataSet_home]
Tpoints = AllPoints
ans_home = []
i=0
random.shuffle(Tpoints)
print(len(Tpoints))
start = time.time()
centerPoint,ans_points = get_max_circle(Tpoints[:300000],r,name="home")
print("time:{}".format(time.time()-start))
print(len(ans_points))
ans_home = []
## 将结果 解析出来
pbar = tqdm(total=5)
for i in range(len(ans_points)):
    p = ans_points[i][0]
    flag = 0
    for j in ans_home:
        if dis(point(j[1],j[0]),p)<0.04:
            flag = 1
            break
    if flag ==1:
        continue
    else:
        count = 0
        for p1 in Tpoints:
            if dis(p1,p)<=0.02:
                count+=1
        ans_home.append([p.y,p.x,count])
        pbar.update(1)
        if len(ans_home)==5:
            break
ans_home.sort(key=lambda x:x[-1],reverse=True)
print(ans_home)


### 测试结果是否正确
for i in range(len(ans_company)):
    for j in range(len(ans_company)):
        if i==j:
            continue
        d = dis(point(ans_company[i][0],ans_company[i][1]),point(ans_company[j][0],ans_company[j][1]))
        if d<=0.04:
            print(d)
for i in range(len(ans_home)):
    for j in range(len(ans_home)):
        if i==j:
            continue
        d = dis(point(ans_home[i][0],ans_home[i][1]),point(ans_home[j][0],ans_home[j][1]))
        if d<=0.04:
            print(d)

# 任务四计算结果
data4 = []
for data in ans_home:
    data4.append(data)
for data in ans_company:
    data4.append(data)
print(data4)
def format_data_list(data):
    return '\n'.join([str(item) for item in data])

def format_data_str(data):
    coor_str = ''
    for item in data:      
        coor_str += ','.join([str(coor) for coor in item]) + '\n'
    return coor_str

def format_data(data1, data2, data3, data4):
    split_str = '\n=====\n'
    # 任务一
    output = ','.join([str(u) for u in data1])
    # 任务二
    output += split_str
    output +=format_data_list(data2)
    # 任务三
    output += split_str
    output +=format_data_list(data3)
    # 任务四
    output += split_str
    output +=format_data_str(data4)
    return output
# 将结果保存至submit.txt文件，之后用于提交
def write_to_txt(data1, data2, data3, data4, path='submit.txt'):
    if len(data1) != 4 or len(data2) != 10 or len(data3) != 8 or len(data4) != 10 :
        return False
    try:
        with open(path, 'w+') as f:
            split_str = '\n=====\n'
            f.write(format_data(data1, data2, data3, data4))
        print('生成提交结果文件成功!')
    except Exception as e:
        print('读写文件失败',e)    

write_to_txt(data1, data2, data3, data4)