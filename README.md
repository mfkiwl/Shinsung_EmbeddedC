# 신성 Embedded C code 입니다.
## 코드 수정사항은 다 여기로 반영해주시면 감사하겠습니다.
### 10월16일데이터: imu20201016_v2.txt    10월21일데이터: imu_V2_201021.txt  
  
    
  
원본 C 코드
======
**10.21.C 코드 검증 완료**

임베디드 코드
======
10.22.9:39 - IMU y,z,축 가속도계 자이로스코프 부호 역전되는 부분 반영, InsKf15_updatelatlonvnve누락 수정  
10.22.12:47 - 필터 조건문수정, 부호역전되는 부분 데이터 취득과정에서 반영, 임베디드로 옮기기 전 원본 C 코드 올림.  
