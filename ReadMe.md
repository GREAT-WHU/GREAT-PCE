# **GREAT-PCE: 武汉大学GREAT团队精密卫星钟差估计软件**

## **软件介绍**

GREAT (GNSS+ REsearch, Application and Teaching) 软件由武汉大学测绘学院设计开发，是一个用于空间大地测量数据处理、精密定位和定轨估钟以及多源融合导航的综合性软件平台。软件中，核心计算模块使用 C++语言(C++17)编写，辅助脚本模块使用 Python3 和 C-Shell 语言实现自动化数据处理。其中，所有 C++模块的编写都遵循 Google 开源项目代码风格指南，并且使用 GIT 工具进行版本控制。GREAT 软件使用 CMake 工具进行编译管理，用户可以灵活选择 GCC、 Clang、MSVC 等主流 C++编译器。目前软件提供了 Windows 和 Linux 平台的命令行应用程序。

GREAT-PCE 是 GREAT 软件中的一个重要模块，主要用于精密卫星钟差解算。GREAT-PCE 由 2 个可移植程序库组成，分别是 LibGREAT 和 LibGnut。LibGREAT 库主要用于最小二乘解算，包括估计中涉及的数据解码、存储以及 PCE 算法的实现。LibGnut 库来源于开源 GNSS 软件 G-nut，包括 GNSS 数据的解码和存储以及基本参数配置模块。

## **软件功能**

本次开源的GREAT-PCE Beta版本主要包括精密卫星钟差估计部分，主要包含以下功能：

- 支持 GPS、GLONASS、Galileo、BDS-2/3 等卫星导航系统

- 支持双频无电离层组合观测值组合方式

- 支持整体解算与仿实时解算两种模式

## **软件包目录结构**

```text
GREAT-PCE_1.0
  ./src	                 源代码 *
    ./app                  GREAT-PCE主程序 *
    ./LibGREAT             精密卫星钟差估计核心算法库 *
    ./LibGnut              Gnut库 *
    ./third-party          第三方库 *
  ./sample_data          算例数据 *
    ./ PCE_2020100         精密卫星钟差估计算例 *
  ./doc                  文档 *
    ./GREAT-PCE_1.0.pdf    GREAT-PCE用户指南 *
```

## **安装和使用**

参见GREAT-PCE _1.0.pdf

## **修改记录**

暂无。

## **参与贡献**

开发人员：

武汉大学GREAT团队, Wuhan University.

三方库：

- GREAT-PCE使用G-Nut库([http://www.pecny.cz](http://www.pecny.cz/)) Copyright (C) 2011-2016 GOP - Geodetic     Observatory Pecny, RIGTC.
- GREAT-PCE使用pugixml库([http://pugixml.org](http://pugixml.org/)) Copyright (C) 2006-2014     Arseny Kapoulkine.
- GREAT-PCE使用Newmat库(http://www.robertnz.net/nm_intro.htm) Copyright (C) 2008: R B     Davies.
- GREAT-MSF使用Eigen库([https://eigen.tuxfamily.org](https://eigen.tuxfamily.org/)) Copyright (C) 2008-2011     Gael Guennebaud

## **致谢**

GREAT-PCE编写过程中部分借鉴参考了葛茂荣老师提供的PANDA软件，再此特别致谢。

## **下载地址**

### GREAT-PCE

GitHub：https://github.com/GREAT-WHU/GREAT-PCE

### 其他GREAT开源软件模块：

GREAT-LAG

- GitHub：https://github.com/GREAT-WHU/GREAT-LAG

GREAT-MSF

- GitHub：https://github.com/GREAT-WHU/GREAT-MSF

GREAT-PVT

- GitHub：https://github.com/GREAT-WHU/GREAT-PVT

欢迎加入QQ群(**1009827379**)参与讨论与交流。

