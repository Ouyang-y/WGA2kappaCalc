# WGA2kappaCalc
从不同长度阵列演化端面图片倒推计算阵列的耦合系数
## 文件结构
``` matlab
WGA2kappaCalc
│
│  WGA2kappaCalc.m  % 主程序
│
├─images
│      Result_init.zip  % init版本结果
│      WGA-13.png       % init版本识别效果
│      Result_V0.1.zip  % V0.1版本结果
│      2_WGA-13.png     % V0.1版本识别效果
├─resource  % 会用到的函数
│      functionSignatures.json  % 自定义代码建议和自动填充
│      gauss2fit.m
│      OSI_rainbow.mat
│      UWGA_evaluate.m
└─test  % 测试文件,说明详见<WGA2kappaCalc.m>
        WGA-1.tiff
        WGA-2.tiff
        WGA-3.tiff
        WGA-4.tiff
        WGA-5.tiff
        WGA-6.tiff
        WGA-7.tiff
        WGA-8.tiff
        WGA-9.tiff
        WGA-10.tiff
        WGA-11.tiff
        WGA-12.tiff
        WGA-13.tiff
        WGA-14.tiff
        wga_01.bmp
        wga_02.bmp
```
## 更新说明
### V0.1
1. 改进注释说明
2. 优化流程逻辑结构
3. 优化拟合光斑方法
    > ![init 识别效果](images/WGA-13.png "init 识别效果")
    > init部分光斑识别不到

    > ![V0.1 识别效果](images/2_WGA-13.png "V0.1 识别效果")
    > V0.1实现改进
4. 优化画图计算逻辑
5. 函数优化并自签名