classdef CRBcal
    properties
        config      % contains the config
        risElementLoc %contias the location of ris elements
        IV    %stores some initial calculations for channel parameters
        FIMch % fisher information matrix for channel parameters
        FIMpo %fisher information matrix for positional parametrs
        PEB %peb vector
        Jacob %jacobian matrix
    end
    
    methods
        function obj = CRBcal(config)
            % constructor: input the config file and the Ris phase profile
            obj.config = config;
            obj.risElementLoc = obj.getRisElLoc;
            obj = obj.getInits;
        end
        
        function RisElLoc = getRisElLoc(obj)
            % calculates a 3XM matrix containing ris element locations
            Lris = (obj.config.Mc-1)*obj.config.risElementDist;
            risXvec = linspace(0,Lris,obj.config.Mc)-Lris/2;
            risZvec = risXvec;
            [risX, risZ] = meshgrid(risXvec,risZvec);
            RisElLoc= [risX(:),zeros(length(risX(:)),1),risZ(:)];
        end
        
        function obj = getInits(obj)
            % calculates some initial values like channel parameters and
            % Power, noise variance , and poistions
            
            obj.IV.sqEs = sqrt(obj.config.TxPower/obj.config.Nsc); %quare root symbol energy
            obj.IV.sqN = sqrt(obj.config.Noise_Factor*obj.config.NPSD*obj.config.Df/2);%square root noise power
            obj.IV.snr = obj.config.TxPower/(obj.config.Nsc*obj.config.Noise_Factor*obj.config.NPSD*obj.config.Df);
            
            obj.IV.uePos = SupPoint(obj.config.xyz(1,:),obj.config.xyz(2,:),obj.config.xyz(3,:)); %ue positions
            obj.IV.BsPos = SupPoint(obj.config.bsPos); %Bs positon
            obj.IV.BsUeVec = SupPoint(obj.IV.uePos.X-obj.IV.BsPos.X, obj.IV.uePos.Y-obj.IV.BsPos.Y, obj.IV.uePos.Z-obj.IV.BsPos.Z);% vector from Bs to Ue
            
            obj.IV.gb = sqrt((obj.config.lambda./(4*pi*obj.IV.BsUeVec.abs(:))).^2).*exp(2*pi*1j*rand(size(obj.IV.BsUeVec.abs(:)))); % los gain
            obj.IV.gr = sqrt((obj.config.lambda./(4*pi*obj.IV.BsPos.abs(:))).^2 .* (obj.config.lambda./(4*pi*obj.IV.uePos.abs(:))).^2).*exp(2*pi*1j*rand(size(obj.IV.BsUeVec.abs(:)))); % reflected path gatin
            obj.IV.tau_b = obj.IV.BsUeVec.abs(:)./obj.config.c; % los ToA
            obj.IV.tau_r = (obj.IV.BsPos.abs(:)+obj.IV.uePos.abs(:))./obj.config.c; % NLOS TOA
            obj.IV.k = obj.IV.uePos.getKv(obj.config.lambda); % wavenumber vector
            obj.IV.kaz = obj.IV.uePos.getKDazv(obj.config.lambda);  %d k /d \phi_az
            obj.IV.kel = obj.IV.uePos.getKDelv(obj.config.lambda);
            obj.IV.RisPhaseProfile = exp(1j*2*pi*rand(obj.config.Mc^2,obj.config.T));
            obj.IV.ur = exp(-1j*obj.IV.k*obj.risElementLoc.') * obj.IV.RisPhaseProfile; %  u_r
            obj.IV.urAz = (exp(-1j*obj.IV.k*obj.risElementLoc.').* (-1j*obj.IV.kaz*obj.risElementLoc.')) * obj.IV.RisPhaseProfile;
            obj.IV.urEl = (exp(-1j*obj.IV.k*obj.risElementLoc.').* (-1j*obj.IV.kel*obj.risElementLoc.')) * obj.IV.RisPhaseProfile;
            
        end
        
        function Rate = getRate(obj)
            % calculates the Spectral Efficiency for random and directional
            % ris phase profiles
            rng(0)
            RisPhaseProfile = exp(1j*2*pi*rand(obj.config.Mc^2,1));
            ur = exp(-1j*obj.IV.k*obj.risElementLoc.') * RisPhaseProfile;
            cyclicPrifix=1.07;
            Noise = (obj.config.Nsc*obj.config.Noise_Factor*obj.config.NPSD*obj.config.Df);
            
            Rate.snr_rand = (obj.config.TxPower*abs(obj.IV.gr.* ur).^2)./Noise;
            Rate.rate_rand = log2(1+Rate.snr_rand)/cyclicPrifix;
            
            Rate.snr_dir = zeros(size(Rate.snr_rand));
            Rate.rate_dir = zeros(size(Rate.rate_rand));
            
            Nn=100;
            for p = 1:size(ur,1)
                P=obj.IV.uePos.getV(1,p);
                peb = obj.PEB(p);
                sig = peb/sqrt(3);
                for n=1:Nn
                    Pn = P + randn(1,3)*sig/sqrt(3);
                    uePos = SupPoint(Pn(1),Pn(2),Pn(3));
                    k  = uePos.getKv(obj.config.lambda);
                    RisPhaseProfile = exp(1j*k*obj.risElementLoc.').';
                    ur = exp(-1j*obj.IV.k(p,:)*obj.risElementLoc.') * RisPhaseProfile;
                    Rate.snr_dir(p) = (obj.config.TxPower*abs(obj.IV.gr(p).*ur).^2)./Noise;
                    Rate.rate_dir(p) = Rate.rate_dir(p)+log2(1+Rate.snr_dir(p))/cyclicPrifix;
                end
                Rate.rate_dir(p)=Rate.rate_dir(p)/Nn;
            end
            
            
        end
        
        
        
        function Jc = getJacobian(obj)
            % returens a cell of the length of number of points each with a 8X8 jacobian matrix
            J = cell(8,8);
            J{1,1} = 1/obj.config.c*(obj.IV.uePos.X-obj.IV.BsPos.X)./(obj.IV.BsUeVec.abs); %d\tau_{b}/d\p1
            J{1,2} = 1/obj.config.c*(obj.IV.uePos.Y-obj.IV.BsPos.Y)./(obj.IV.BsUeVec.abs); %d\tau_{b}/d\p2
            J{1,3} = 1/obj.config.c*(obj.IV.uePos.Z-obj.IV.BsPos.Z)./(obj.IV.BsUeVec.abs); %d\tau_{b}/d\p3
            J{2,1} = 1/obj.config.c*(obj.IV.uePos.X)./(obj.IV.uePos.abs); %d\tau_{r}/d\p1
            J{2,2} = 1/obj.config.c*(obj.IV.uePos.Y)./(obj.IV.uePos.abs); %d\tau_{r}/d\p2
            J{2,3} = 1/obj.config.c*(obj.IV.uePos.Z)./(obj.IV.uePos.abs); %d\tau_{r}/d\p3
            J{3,1} = -obj.IV.uePos.Y./(obj.IV.uePos.XY2); %d\phi_{az}/d\p1
            J{3,2} = obj.IV.uePos.X./(obj.IV.uePos.XY2); %d\phi_{az}/d\p1
            J{4,1} = obj.IV.uePos.X.*obj.IV.uePos.Z./(obj.IV.uePos.abs2.*obj.IV.uePos.XY); %d\phi_{el}/d\p1
            J{4,2} = obj.IV.uePos.Y.*obj.IV.uePos.Z./(obj.IV.uePos.abs2.*obj.IV.uePos.XY); %d\phi_{el}/d\p2
            J{4,3} = (-obj.IV.uePos.abs2 + obj.IV.uePos.Z.*obj.IV.uePos.Z)./(obj.IV.uePos.abs2.*obj.IV.uePos.XY); %d\phi_{el}/d\p3
            Jc = cell(size(obj.IV.uePos.X));
            for ir=1:size(obj.IV.uePos.X,1)
                for ic=1:size(obj.IV.uePos.X,2)
                    Jc{ir,ic}=[J{1,1}(ir,ic),J{1,2}(ir,ic),J{1,3}(ir,ic),1,0,0,0,0;...
                        J{2,1}(ir,ic),J{2,2}(ir,ic),J{2,3}(ir,ic),1,0,0,0,0;...
                        J{3,1}(ir,ic),J{3,2}(ir,ic),0,0,0,0,0,0;...
                        J{4,1}(ir,ic),J{4,2}(ir,ic),J{4,3}(ir,ic),0,0,0,0,0;...
                        0,0,0,0,1,0,0,0;0,0,0,0,0,1,0,0;0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,1];
                end
            end
        end
        
        function FIMch = getFIMCh(obj)
            % calculates the fisher information matrix for channel
            % parameters
            db = exp(-1j*2*pi*mod(obj.config.Df*obj.IV.tau_b*[0:obj.config.Nsc-1],1));
            dr = exp(-1j*2*pi*mod(obj.config.Df*obj.IV.tau_r*[0:obj.config.Nsc-1],1));
            d = -1j*2*pi*[0:obj.config.Nsc-1]*obj.config.Df;
            FIMch = cell(1,size(obj.config.xyz,2));
            tic
            for ip = 1:size(obj.config.xyz,2)
                fprintf('Calculating CRB %.2f %% done! \n', (ip-1)/size(obj.config.xyz,2)*100);
                FIMch{1,ip} = zeros(8);
                for it = 1:obj.config.T
                    for in = 1:obj.config.Nsc
                        u = [
                            d(in)*obj.IV.gb(ip)*db(ip,in),...
                            d(in)*obj.IV.gr(ip)*dr(ip,in)*obj.IV.ur(ip,it),...
                            obj.IV.gr(ip)*dr(ip,in)*obj.IV.urAz(ip,it),...
                            obj.IV.gr(ip)*dr(ip,in)*obj.IV.urEl(ip,it),...
                            db(ip,in),...
                            1j*db(ip,in),...
                            dr(ip,in)*obj.IV.ur(ip,it),...
                            1j*dr(ip,in)*obj.IV.ur(ip,it),...
                            ];
                        FIMch{1,ip} = FIMch{1,ip} + 2*obj.IV.snr* real(u'*u);
                    end
                end
                ElTime = toc;
                fprintf('Est. remaining time:  %.2f minuites \n', ElTime/ip*(size(obj.config.xyz,2)-ip)/60);
                
            end
        end
        
        function obj = PEBcalc(obj)
            %calculates fisher information for positional parameters and
            %PEB
            obj.Jacob = obj.getJacobian;
            obj.FIMch = obj.getFIMCh;
            obj.FIMpo = cell(1,size(obj.config.xyz,2));
            obj.PEB = zeros(1,size(obj.config.xyz,2));
            for ip = 1:size(obj.config.xyz,2)
                Jcb = obj.Jacob{1,ip};
                obj.FIMpo{1,ip} = Jcb.' * obj.FIMch{1,ip} * Jcb;
                
                fim = obj.FIMpo{1,ip};% using Schur’s complement
                A = fim(1:4,1:4);
                B = fim(1:4,5:8);
                C = fim(5:8,1:4);
                D = fim(5:8,5:8);
                eqFim = A-B*inv(D)*C;
                FimInv = inv(eqFim);
                obj.PEB(1,ip) = sqrt(trace(FimInv(1:3,1:3)));
            end
        end
        
        
        
    end
end


