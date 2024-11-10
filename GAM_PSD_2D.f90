    Program PSD_2D
    Implicit None
    Integer,Allocatable::Image(:,:)
    Real,Allocatable::Radius_2D(:,:)
    Real,Allocatable::Radius_2D_5(:,:)
    Real,Allocatable::Radius_2D_0(:,:)
    Real,Allocatable::Radius_Max_5(:,:)
    Real,Allocatable::Radius_Max_0(:,:)
    Real,Allocatable::Dist_5(:,:)
    Real,Allocatable::Dist_0(:,:)
    Real,Allocatable::PSD(:,:)
    Integer Height,Width,Radius,Num_Pore,Num_Radius
    Integer i,j,r
    Integer Time(8)
    Open(1,File='Parameters.dat')
    Read(1,*) Height,Width,Radius
    Open(2,File='Image.dat')
    Allocate(Image(Height,Width))
    Read(2,*) ((Image(i,j),j=1,Width),i=1,Height)
    Allocate(Dist_5(2*Radius-1,2*Radius-1))
    Allocate(Dist_0(2*Radius,2*Radius))
    !$OMP PARALLEL
    !$OMP DO
    Do j=1,2*Radius-1
        Do i=1,2*Radius-1
            Dist_5(i,j)=((i-Radius)**2+(j-Radius)**2)**0.5
        End Do
    End Do
    !$OMP END DO
    !$OMP END PARALLEL
    !$OMP PARALLEL
    !$OMP DO
    Do j=1,2*Radius
        Do i=1,2*Radius
            Dist_0(i,j)=((i-Radius-0.5)**2+(j-Radius-0.5)**2)**0.5
        End Do
    End Do
    !$OMP END DO
    !$OMP END PARALLEL
    Allocate(Radius_2D(Height,Width))
    Radius_2D=0
    Allocate(Radius_2D_5(Height,Width))
    Radius_2D_5=0
    Allocate(Radius_2D_0(Height,Width))
    Radius_2D_0=0
    Allocate(Radius_Max_5(Height,Width))
    Radius_Max_5=0
    Where(Image==1) Radius_Max_5=0.5
    Allocate(Radius_Max_0(Height,Width))
    Radius_Max_0=0
    Open(3,File='Date_And_Time.dat')
    Call Date_And_Time(Values=Time)
    Write(3,20) Time(1),"/",Time(2),"/",Time(3)," ",Time(5),":",Time(6),":",Time(7)
    !$OMP PARALLEL
    !$OMP DO
    Do j=1,Width
        Do i=1,Height
            If (Image(i,j)==1) Then
                DO r=2,Radius
                    If ((i-r+1>=1).AND.(i+r-1<=Height).AND.(j-r+1>=1).AND.(j+r-1<=Width)) Then
                        If (Minval(Image(i-r+1:i+r-1,j-r+1:j+r-1)+(Dist_5(Radius-r+1:Radius+r-1,Radius-r+1:Radius+r-1)<=(r-0.5)))==0) Then
                            Radius_Max_5(i,j)=r-0.5
                        Else
                            Exit
                        End If
                    Else
                        Exit
                    End If
                End Do
                DO r=1,Radius
                    If ((i-r+1>=1).AND.(i+r<=Height).AND.(j-r+1>=1).AND.(j+r<=Width)) Then
                        If (Minval(Image(i-r+1:i+r,j-r+1:j+r)+(Dist_0(Radius-r+1:Radius+r,Radius-r+1:Radius+r)<=r))==0) Then
                            Radius_Max_0(i,j)=r
                        Else
                            Exit
                        End If
                    Else
                        Exit
                    End If
                End Do
            End If
        End Do
    End Do
    !$OMP END DO
    !$OMP END PARALLEL
    !$OMP PARALLEL
    !$OMP DO
    Do j=1,Width
        Do i=1,Height
            If (Radius_Max_5(i,j)/=0) Then
                Where(((Dist_5(Radius-Ceiling(Radius_Max_5(i,j))+1:Radius+Ceiling(Radius_Max_5(i,j))-1,Radius-Ceiling(Radius_Max_5(i,j))+1:Radius+Ceiling(Radius_Max_5(i,j))-1)<=Radius_Max_5(i,j))==-1).AND.(Radius_2D_5(i-Ceiling(Radius_Max_5(i,j))+1:i+Ceiling(Radius_Max_5(i,j))-1,j-Ceiling(Radius_Max_5(i,j))+1:j+Ceiling(Radius_Max_5(i,j))-1)<Radius_Max_5(i,j)))
                    Radius_2D_5(i-Ceiling(Radius_Max_5(i,j))+1:i+Ceiling(Radius_Max_5(i,j))-1,j-Ceiling(Radius_Max_5(i,j))+1:j+Ceiling(Radius_Max_5(i,j))-1)=Radius_Max_5(i,j)
                End Where
            End If
            If (Radius_Max_0(i,j)/=0) Then
                Where(((Dist_0(Radius-Ceiling(Radius_Max_0(i,j))+1:Radius+Ceiling(Radius_Max_0(i,j)),Radius-Ceiling(Radius_Max_0(i,j))+1:Radius+Ceiling(Radius_Max_0(i,j)))<=Radius_Max_0(i,j))==-1).AND.(Radius_2D_0(i-Ceiling(Radius_Max_0(i,j))+1:i+Ceiling(Radius_Max_0(i,j)),j-Ceiling(Radius_Max_0(i,j))+1:j+Ceiling(Radius_Max_0(i,j)))<Radius_Max_0(i,j)))
                    Radius_2D_0(i-Ceiling(Radius_Max_0(i,j))+1:i+Ceiling(Radius_Max_0(i,j)),j-Ceiling(Radius_Max_0(i,j))+1:j+Ceiling(Radius_Max_0(i,j)))=Radius_Max_0(i,j)
                End Where
            End If
        End Do
    End Do
    !$OMP END DO
    !$OMP END PARALLEL
    !$OMP PARALLEL
    !$OMP DO
    Do j=1,Width
        Do i=1,Height
            If (Radius_2D_0(i,j)>Radius_2D_5(i,j)) Then
                Radius_2D(i,j)=Radius_2D_0(i,j)
            Else
                Radius_2D(i,j)=Radius_2D_5(i,j)
            End If
        End Do
    End Do
    !$OMP END DO
    !$OMP END PARALLEL
    Call Date_And_Time(Values=Time)
    Write(3,20) Time(1),"/",Time(2),"/",Time(3)," ",Time(5),":",Time(6),":",Time(7)
    Num_Pore=Sum(Image)
    Num_Radius=Int(Maxval(Radius_2D)*2)
    Allocate(PSD(Num_Radius,4))
    Do i=1,Num_Radius
        PSD(i,1)=i*0.5
        PSD(i,2)=Count(Radius_2D==i*0.5)
        PSD(i,3)=(PSD(i,2)*1.0)/(Num_Pore*1.0)*100.0
        PSD(i,4)=Sum(PSD(1:i,3))
    End Do
    Open(4,File='PSD.dat')
    Write(4,10) ((PSD(i,j),j=1,4),i=1,Num_Radius)
10  Format(F5.1,F10.0,F10.4,F10.4)
20  Format(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)
    End Program PSD_2D