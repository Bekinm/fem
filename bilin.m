function [fs_next, ue_next, Ehys_next, kT_next] = bilin(k1,k2,fy,u,ue,fs,du,Ehys)
        fs_next=fs+k1*du;
        u_next=u+du;
        if fs_next>k2*u_next+fy*(1-k2/k1)
            fs_next=k2*u_next+fy*(1-k2/k1);
            ue_next=ue+((k2*u+fy*(1-k2/k1)-fs)/(k1-k2));
            Ehys_next=Ehys+abs((fs_next+fs)*0.5*(du-((k2*u+fy*(1-k2/k1)-fs)/(k1-k2))));
            kT_next=k2;
        elseif fs_next<k2*u_next-fy*(1-k2/k1)
            fs_next=k2*u_next-fy*(1-k2/k1);
            ue_next=ue+((k2*u-fy*(1-k2/k1)-fs)/(k1-k2));
            Ehys_next=Ehys+abs((fs_next+fs)*0.5*(du-((k2*u-fy*(1-k2/k1)-fs)/(k1-k2))));
            kT_next=k2;
        else
            ue_next=ue+du;
            Ehys_next=Ehys;
            kT_next=k1;
        end
end