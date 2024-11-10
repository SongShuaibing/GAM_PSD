    Program PSD_3D
    Implicit None
    Integer,Allocatable::Image(:,:,:)
    Real,Allocatable::Radius_3D(:,:,:)
    Real,Allocatable::Radius_3D_5(:,:,:)
    Real,Allocatable::Radius_3D_0(:,:,:)
    Real,Allocatable::Radius_Max_5(:,:,:)
    Real,Allocatable::Radius_Max_0(:,:,:)
    Real,Allocatable::Dist_5(:,:,:)
    Real,Allocatable::Dist_0(:,:,:)
    Real,Allocatable::PSD(:,:)
    Integer Height,Width,Depth,Radius,Num_Pore,Num_Radius
    Integer i,j,k,r
    Integer Time(8)
    Open(1,File='Parameters.dat')
    Read(1,*) Height,Width,Depth,Radius
    Open(2,File='Image.dat')
    Allocate(Image(Height,Width,Depth))
    Read(2,*) (((Image(i,j,k),j=1,Width),i=1,Height),k=1,Depth)
    Allocate(Dist_5(2*Radius-1,2*Radius-1,2*Radius-1))
    Allocate(Dist_0(2*Radius,2*Radius,2*Radius))
    !$OMP PARALLEL
    !$OMP DO
    Do k=1,2*Radius-1
        Do j=1,2*Radius-1
            Do i=1,2*Radius-1
                Dist_5(i,j,k)=((i-Radius)**2+(j-Radius)**2+(k-Radius)**2)**0.5
            End Do
        End Do
    End Do
    !$OMP END DO
    !$OMP END PARALLEL
    !$OMP PARALLEL
    !$OMP DO
    Do k=1,2*Radius
        Do j=1,2*Radius
            Do i=1,2*Radius
                Dist_0(i,j,k)=((i-Radius-0.5)**2+(j-Radius-0.5)**2+(k-Radius-0.5)**2)**0.5
            End Do
        End Do
    End Do
    !$OMP END DO
    !$OMP END PARALLEL
    Allocate(Radius_3D(Height,Width,Depth))
    Radius_3D=0
    Allocate(Radius_3D_5(Height,Width,Depth))
    Radius_3D_5=0
    Allocate(Radius_3D_0(Height,Width,Depth))
    Radius_3D_0=0
    Allocate(Radius_Max_5(Height,Width,Depth))
    Radius_Max_5=0
    Where(Image==1) Radius_Max_5=0.5
    Allocate(Radius_Max_0(Height,Width,Depth))
    Radius_Max_0=0
    Open(3,File='Date_And_Time.dat')
    Call Date_And_Time(Values=Time)
    Write(3,20) Time(1),"/",Time(2),"/",Time(3)," ",Time(5),":",Time(6),":",Time(7)
    !$OMP PARALLEL
    !$OMP DO
    Do k=1,Depth
        Do j=1,Width
            Do i=1,Height
                If (Image(i,j,k)==1) Then
                    DO r=2,Radius
                        If ((i-r+1>=1).AND.(i+r-1<=Height).AND.(j-r+1>=1).AND.(j+r-1<=Width).AND.(k-r+1>=1).AND.(k+r-1<=Depth)) Then
                            If (Minval(Image(i-r+1:i+r-1,j-r+1:j+r-1,k-r+1:k+r-1)+(Dist_5(Radius-r+1:Radius+r-1,Radius-r+1:Radius+r-1,Radius-r+1:Radius+r-1)<=(r-0.5)))==0) Then
                                Radius_Max_5(i,j,k)=r-0.5
                            Else
                                Exit
                            End If
                        Else
                            Exit
                        End If
                    End Do
                    DO r=1,Radius
                        If ((i-r+1>=1).AND.(i+r<=Height).AND.(j-r+1>=1).AND.(j+r<=Width).AND.(k-r+1>=1).AND.(k+r<=Depth)) Then
                            If (Minval(Image(i-r+1:i+r,j-r+1:j+r,k-r+1:k+r)+(Dist_0(Radius-r+1:Radius+r,Radius-r+1:Radius+r,Radius-r+1:Radius+r)<=r))==0) Then
                                Radius_Max_0(i,j,k)=r
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
    End Do
    !$OMP END DO
    !$OMP END PARALLEL
    !$OMP PARALLEL
    !$OMP DO
    Do k=1,Depth
        Do j=1,Width
            Do i=1,Height
                If (Radius_Max_5(i,j,k)/=0) Then
                    Where(((Dist_5(Radius-Ceiling(Radius_Max_5(i,j,k))+1:Radius+Ceiling(Radius_Max_5(i,j,k))-1,Radius-Ceiling(Radius_Max_5(i,j,k))+1:Radius+Ceiling(Radius_Max_5(i,j,k))-1,Radius-Ceiling(Radius_Max_5(i,j,k))+1:Radius+Ceiling(Radius_Max_5(i,j,k))-1)<=Radius_Max_5(i,j,k))==-1).AND.(Radius_3D_5(i-Ceiling(Radius_Max_5(i,j,k))+1:i+Ceiling(Radius_Max_5(i,j,k))-1,j-Ceiling(Radius_Max_5(i,j,k))+1:j+Ceiling(Radius_Max_5(i,j,k))-1,k-Ceiling(Radius_Max_5(i,j,k))+1:k+Ceiling(Radius_Max_5(i,j,k))-1)<Radius_Max_5(i,j,k)))
                        Radius_3D_5(i-Ceiling(Radius_Max_5(i,j,k))+1:i+Ceiling(Radius_Max_5(i,j,k))-1,j-Ceiling(Radius_Max_5(i,j,k))+1:j+Ceiling(Radius_Max_5(i,j,k))-1,k-Ceiling(Radius_Max_5(i,j,k))+1:k+Ceiling(Radius_Max_5(i,j,k))-1)=Radius_Max_5(i,j,k)
                    End Where
                End If
                If (Radius_Max_0(i,j,k)/=0) Then
                    Where(((Dist_0(Radius-Ceiling(Radius_Max_0(i,j,k))+1:Radius+Ceiling(Radius_Max_0(i,j,k)),Radius-Ceiling(Radius_Max_0(i,j,k))+1:Radius+Ceiling(Radius_Max_0(i,j,k)),Radius-Ceiling(Radius_Max_0(i,j,k))+1:Radius+Ceiling(Radius_Max_0(i,j,k)))<=Radius_Max_0(i,j,k))==-1).AND.(Radius_3D_0(i-Ceiling(Radius_Max_0(i,j,k))+1:i+Ceiling(Radius_Max_0(i,j,k)),j-Ceiling(Radius_Max_0(i,j,k))+1:j+Ceiling(Radius_Max_0(i,j,k)),k-Ceiling(Radius_Max_0(i,j,k))+1:k+Ceiling(Radius_Max_0(i,j,k)))<Radius_Max_0(i,j,k)))
                        Radius_3D_0(i-Ceiling(Radius_Max_0(i,j,k))+1:i+Ceiling(Radius_Max_0(i,j,k)),j-Ceiling(Radius_Max_0(i,j,k))+1:j+Ceiling(Radius_Max_0(i,j,k)),k-Ceiling(Radius_Max_0(i,j,k))+1:k+Ceiling(Radius_Max_0(i,j,k)))=Radius_Max_0(i,j,k)
                    End Where
                End If
            End Do
        End Do
    End Do
    !$OMP END DO
    !$OMP END PARALLEL
    !$OMP PARALLEL
    !$OMP DO
    Do k=1,Depth
        Do j=1,Width
            Do i=1,Height
                If (Radius_3D_0(i,j,k)>Radius_3D_5(i,j,k)) Then
                    Radius_3D(i,j,k)=Radius_3D_0(i,j,k)
                Else
                    Radius_3D(i,j,k)=Radius_3D_5(i,j,k)
                End If
            End Do
        End Do
    End Do
    !$OMP END DO
    !$OMP END PARALLEL
    Call Date_And_Time(Values=Time)
    Write(3,20) Time(1),"/",Time(2),"/",Time(3)," ",Time(5),":",Time(6),":",Time(7)
    Num_Pore=Sum(Image)
    Num_Radius=Int(Maxval(Radius_3D)*2)
    Allocate(PSD(Num_Radius,4))
    Do i=1,Num_Radius
        PSD(i,1)=i*0.5
        PSD(i,2)=Count(Radius_3D==i*0.5)
        PSD(i,3)=(PSD(i,2)*1.0)/(Num_Pore*1.0)*100.0
        PSD(i,4)=Sum(PSD(1:i,3))
    End Do
    Open(4,File='PSD.dat')
    Write(4,10) ((PSD(i,j),j=1,4),i=1,Num_Radius)
10  Format(F5.1,F10.0,F10.4,F10.4)
20  Format(I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)
    End Program PSD_3D